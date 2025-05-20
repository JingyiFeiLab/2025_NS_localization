% Initialize a new reader per worker as Bio-Formats is not thread safe
reader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
reader.setId([rawimage_folder filesep wellname '.nd2']);
seriesCount = reader.getSeriesCount();

cell_raw=cell(100000,4);
cell_raw_ins=zeros(100000,17);
bg_median=zeros(seriesCount,2);
% cell_raw (cell): cells near edges removed, without quality control, bg included
% 1 : I_raw_crop
% 2: bw_cell_crop
% 3: bw_nu
% 4: cell_id_all
% cell_raw_ins (mat) : 
% 1: cell_id_all
% 2: # of series
% 3: maskID in images
% 4: IR_cell
% 5: IG_cell
% 6: IR_nu
% 7: IG_nu
% 8: IR_cyto
% 9: IG_cyto
% 10: A_cell
% 11: A_nu
% 12: cR_cell
% 13: cG_cell
% 14: cR_nu
% 15: cG_nu
% 16: cR_cyto
% 17: cG_cyto
% bg_info: 1: bg_R 2:bg_G

cell_id_all=0;
for i=1:seriesCount
    reader.setSeries(i - 1); % note setSeries starts from 0
    %iPlane = reader.getIndex(iZ-1, iC -1, iT-1) + 1;
    % note getIndex starts from 0 but --- iPlane starts from 1
    I_raw = uint16(zeros(1024,1024,3));
    I_raw(:,:,1) = bfGetPlane(reader, 1);
    I_raw(:,:,2) = bfGetPlane(reader, 2);
    I_raw(:,:,3) = bfGetPlane(reader, 3);

    pathTomask=[outline_folder filesep 'ad-' wellname '_' num2str(i) '_cp_masks.png'];
    if isfile(pathTomask) == 0
        continue
    end
    mask=imread(pathTomask);

    cell_num_image=max(max(mask));
    bw_bg=(mask==0);
    %Note now bg is estimated by medians of intensities
    bg_pixels_R=nonzeros(I_raw(:,:,1));
    bg_pixels_G=nonzeros(I_raw(:,:,2));
    bg_median(i,1)=median(bg_pixels_R);
    bg_median(i,2)=median(bg_pixels_G);

    for j=1:cell_num_image
        bw_cell=(mask==j);
        [rows,cols]=find(bw_cell);
        if isempty(rows)
            continue
        end
        mincol=min(cols);
        maxcol=max(cols);
        minrow=min(rows);
        maxrow=max(rows);

        if mincol==1 || maxcol==size(mask,1) || minrow==1 || maxrow==size(mask,1)
            continue
        end

        cell_id_all=cell_id_all+1;
        bw_cell_crop=bw_cell(minrow:maxrow,mincol:maxcol);
        I_crop=I_raw(minrow:maxrow,mincol:maxcol,:);

        cell_raw{cell_id_all,1}=I_crop;
        cell_raw{cell_id_all,2}=bw_cell_crop;
        cell_raw{cell_id_all,4}=cell_id_all;
        cell_raw_ins(cell_id_all,1)=cell_id_all;
        cell_raw_ins(cell_id_all,2)=i;
        cell_raw_ins(cell_id_all,3)=j;

    end
end

reader.close()
cell_raw=cell_raw(~cellfun('isempty',cell_raw(:,1)),:);
bg_median(bg_median(:,1)==0,:)=[];
cell_raw_ins(cell_raw_ins(:,1)==0,:)=[];
bg_info=median(bg_median);
save([save_folder filesep wellcons 'cell_raw.mat'],'cell_raw')
save([save_folder filesep wellcons 'cell_raw_ins.mat'],'cell_raw_ins','bg_info')