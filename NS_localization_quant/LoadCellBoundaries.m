%% LoadCellBoundaries.m
% Features:
% 1) Bio-Formats importer: 
% direct access to .nd2 files , Bio-image information extracted from standardized OME metadata
% 2) 2 modes: manual (by manual ND acquisition), auto (by automatic well-imaging)
% 3) To find the focus plane: Detect the in-focus plane by absolute Laplacian
% Outputs description:
% 1) 'cell_raw' is a cell array with the following columns:
% - col1: 'I_raw_crop' - The cropped raw image data for the cell.
% - col2: 'bw_cell_crop' - Binary mask of the cropped cell (from cellpose).
% - col3: 'imDAPI_crop' - Cropped DAPI image.
% - col4: 'cell_id_image' - The unique identifier for the cell within the image.
% - col5: 'image_id' - The identifier for the image the cell belongs to.
% - col6: 'cell_id_all' - The global identifier for the cell across all images.
% - col7: 'z_max' - The z-index for focus plane.
% - col8: 'cytoplasm_edge' - Indicates if the cell is near the edge (1 for true, 0 for false).
% 2)'cell_intensity' is a matrix with the following columns:
% - col1: Area: The area of the cell.
% - col2: total cellular intensity of the first channel.
% - col3: total cellular intensity of the second channel.
% - col4: total cellular intensity of the third channel.
% 3) back_R, back_G, back_B: camera background (from individual images) estimated
% by regions unoccupied by cells.

%%
% Pre-allocate memory for variables
cell_intensity=zeros(50000,4);
cell_raw=cell(50000,8);
back_R=zeros(2000,1);
back_G=zeros(2000,1);
back_B=zeros(2000,1);

cell_id_all=0;
num_of_image=0;

if strcmp(image_type,'manual')
    FileList = dir([rawimage_folder filesep folder filesep subfolder]);
    maxjj=length(FileList);
else
    maxjj=1;
end

for jj=1:maxjj
    if strcmp(image_type,'manual')
        if FileList(jj).isdir
            continue
        end
        pathToFile=[rawimage_folder filesep folder filesep subfolder filesep FileList(jj).name];
        reader = bfGetReader(pathToFile);
    else
        pathToFile = [rawimage_folder filesep folder filesep subfolder filesep image_name '.nd2'];
        reader = bfGetReader();
        reader = loci.formats.Memoizer(reader);
        reader.setId(pathToFile);
    end
    seriesCount = reader.getSeriesCount();
    zstackCount=reader.getSizeZ();

    for i=1:seriesCount
        num_of_image=num_of_image+1;
        pathTomask=[outline_folder filesep folder filesep subfolder filesep 'ad_' num2str(num_of_image) '_cp_masks.png'];
        if isfile(pathTomask) == 0
            continue
        end
        mask=imread(pathTomask);
        disp(['Loading Masks: Image ' num2str(num_of_image)])
        I_R=uint16(zeros(1024,1024,zstackCount));
        I_G=uint16(zeros(1024,1024,zstackCount));
        I_B=uint16(zeros(1024,1024,zstackCount));
        imDAPI=uint16(zeros(1024,1024,zstackCount));
        reader.setSeries(i - 1);
        for ii=1:zstackCount
            iPlane_R = reader.getIndex(ii-1, R_channel-1, 1-1) + 1;
            iPlane_G = reader.getIndex(ii-1, G_channel-1, 1-1) + 1;
            iPlane_B = reader.getIndex(ii-1, B_channel-1, 1-1) + 1;
            iPlane_DAPI = reader.getIndex(ii-1, DAPI_channel-1, 1-1) + 1;
            I_R(:,:,ii) = bfGetPlane(reader, iPlane_R);
            I_G(:,:,ii) = bfGetPlane(reader, iPlane_G);
            I_B(:,:,ii) = bfGetPlane(reader, iPlane_B);
            imDAPI(:,:,ii) = bfGetPlane(reader, iPlane_DAPI);
        end

        % camera background (from individual images) estimated by regions
        % unoccupied by cells. 
        IR_elements = reshape(uint16(mask==0).*max(I_R, [], 3),[],1);
        IR_elements = nonzeros(IR_elements);
        IG_elements = reshape(uint16(mask==0).*max(I_G, [], 3),[],1);
        IG_elements = nonzeros(IG_elements);
        IB_elements = reshape(uint16(mask==0).*max(I_B, [], 3),[],1);
        IB_elements = nonzeros(IB_elements);
        back_R(num_of_image,1)=median(IR_elements);
        back_G(num_of_image,1)=median(IG_elements);
        back_B(num_of_image,1)=median(IB_elements);

        cell_num_image=max(max(mask));

        for j=1:cell_num_image
            bw_cell=(mask==j);
            [rows,cols]=find(bw_cell);
            if isempty(rows)
                continue
            end

            cell_id_all=cell_id_all+1;
            mincol=min(cols);
            maxcol=max(cols);
            minrow=min(rows);
            maxrow=max(rows);
            cell_raw{cell_id_all,8}=0;
            if mincol==1 || maxcol==size(mask,1) || minrow==1 || maxrow==size(mask,1)
                cell_raw{cell_id_all,8}=1;
            end

            bw_cell_crop=bw_cell(minrow:maxrow,mincol:maxcol);

            % crop the image into multi-stack I_R, I_G, I_B, imDAPI
            I_R_crop=I_R(minrow:maxrow,mincol:maxcol,:).*uint16(bw_cell_crop);
            I_G_crop=I_G(minrow:maxrow,mincol:maxcol,:).*uint16(bw_cell_crop);
            I_B_crop=I_B(minrow:maxrow,mincol:maxcol,:).*uint16(bw_cell_crop);
            imDAPI_crop=imDAPI(minrow:maxrow,mincol:maxcol,:).*uint16(bw_cell_crop);

            % find the focus plane by computing absolute Laplacian
            AL=zeros(zstackCount,1);
            for z=1:zstackCount
                I=imgaussfilt(I_B_crop(:,:,z),0.5);
                % AL: absolute Laplacian
                Ixy=I(2:end-1,2:end-1);
                Ixm1y=I(1:end-2,2:end-1);
                Ixp1y=I(3:end,2:end-1);
                Ixym1=I(2:end-1,1:end-2);
                Ixyp1=I(2:end-1,3:end);
                AL(z)=sum(sum(abs(2*Ixy-Ixm1y-Ixp1y)))+sum(sum(abs(2*Ixy-Ixym1-Ixyp1)));
                AL(z)=AL(z)/(size(I,1)*size(I,2));
            end
            [~,z_max]=max(AL);

            I_raw=zeros(size(I_B_crop,1),size(I_B_crop,2),'uint16');
            I_raw(:,:,1)=I_R_crop(:,:,z_max);
            I_raw(:,:,2)=I_G_crop(:,:,z_max);
            I_raw(:,:,3)=I_B_crop(:,:,z_max);
            imDAPI_crop_max=imDAPI_crop(:,:,z_max);

            % Alternatively, use maximum projection instead of finding the in-focus plane
            % I_raw(:,:,1)=max(I_R_crop, [], 3);
            % I_raw(:,:,2)=max(I_G_crop, [], 3);
            % I_raw(:,:,3)=max(I_B_crop, [], 3);
            % imDAPI_crop_max=max(imDAPI_crop,[],3);

            % Generate statistics
            cell_intensity(cell_id_all,2)=sum(sum(I_raw(:,:,1)));
            cell_intensity(cell_id_all,3)=sum(sum(I_raw(:,:,2)));
            cell_intensity(cell_id_all,4)=sum(sum(I_raw(:,:,3)));
            cell_intensity(cell_id_all,1)=sum(sum(uint16(bw_cell_crop)));
            cell_raw{cell_id_all,1}=I_raw;
            cell_raw{cell_id_all,2}=bw_cell_crop;
            cell_raw{cell_id_all,3}=imDAPI_crop_max;
            cell_raw{cell_id_all,4}=j;
            cell_raw{cell_id_all,5}=num_of_image;
            cell_raw{cell_id_all,6}=cell_id_all;
            cell_raw{cell_id_all,7}=z_max;
        end
    end

end

cell_intensity(cell_intensity(:,2)==0,:)=[];
cell_raw=cell_raw(~cellfun('isempty',cell_raw(:,1)),:);
back_R(back_R==0,:)=[];
back_G(back_G==0,:)=[];
back_B(back_B==0,:)=[];
cam_bg_R=median(double(back_R));
cam_bg_G=median(double(back_G));
cam_bg_B=median(double(back_B));
save([cellname filesep consname 'cell_intensity.mat'],'cell_intensity')
save([cellname filesep consname 'cell_raw.mat'],'cell_raw')
save([cellname filesep consname 'cam_bg.mat'],'cam_bg_G','cam_bg_R','cam_bg_B')
