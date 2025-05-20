load([cellname filesep consname 'cell_raw.mat'])
load([cellname filesep consname 'cell_intensity.mat'])
load([nucleusname filesep consname 'nu_intensity.mat'])
load([nucleusname filesep consname 'nu_raw.mat'])
load([cellname filesep consname 'cam_bg.mat'])
cell_num_all=size(nu_intensity,1);
% Outputs descriptions:
% trans(matrix):
% col1: cell_id_image - Identifier for the cell in the image
% col2: image_id - Identifier for the corresponding image
% col3: cell_id_all - Identifier for the cell across all images
% col4: z_max - The focus plane z-value
% untrans(matrix):
% col1: cell_id_image - Identifier for the cell in the image
% col2: image_id - Identifier for the corresponding image
% col3: cell_id_all - Identifier for the cell across all images
% col4: z_max - The focus plane z-value
% bg(matrix):
% Row 1: [cam_bg_R cam_bg_G cam_bg_B] - camera background (from individual
% images) estimated from LoadCellBoundaries.m
% Row 2: [nu_R nu_G nu_B] - Background color for nuclei in red, green, and blue channels
% Row 3: [cyto_R cyto_G cyto_B] - Background color for cytoplasm in red, green, and blue channels

trans=zeros(cell_num_all,4);
untrans=zeros(cell_num_all,4);
bg=zeros(3,3);
c_nu_R=nu_intensity(:,2)./nu_intensity(:,1);

for i=1:cell_num_all
    % filter out cells without nucleus
    if nu_intensity(i,1)==0
        continue
    end
    % filter out cells with oversaturated pixels
    I_raw=nu_raw{i,1};
    sat_R=length(find(I_raw(:,:,1)>4094));
    sat_G=length(find(I_raw(:,:,2)>4094));
    if sat_R>sat_num_threshold_R || sat_G>sat_num_threshold_G
        continue
    end
    % choose transfected cells

    I_nu_R=nu_raw{i,1}(:,:,1).*uint16(nu_raw{i,2});
    %     if isnan(c_nu_R(i))
    %         continue
    %     end

    % There are two types of input for `trans_threshold`:
    % 1. If trans_threshold < 1: It represents a percentage (e.g., 0.5).
    % 2. If trans_threshold > 1: It represents intensity values (e.g., pixel intensity thresholds).
    if trans_threshold < 1
        transfected_threshold=nanmedian(c_nu_R)+trans_threshold*nanmedian(c_nu_R);
    end
    if trans_threshold > 1
        transfected_threshold=trans_threshold;
    end

    if c_nu_R(i)>transfected_threshold
        trans(i,1)=cell_raw{i,4};
        trans(i,2)=cell_raw{i,5};
        trans(i,3)=cell_raw{i,6};
        trans(i,4)=cell_raw{i,7};
    else
        untrans(i,1)=cell_raw{i,4};
        untrans(i,2)=cell_raw{i,5};
        untrans(i,3)=cell_raw{i,6};
        untrans(i,4)=cell_raw{i,7};
    end
end
trans(trans(:,1)==0,:)=[];
untrans(untrans(:,1)==0,:)=[];
untrans_id=untrans(:,3);

% background estimated from
bg(1,1)=cam_bg_R;
bg(1,2)=cam_bg_R;
bg(1,3)=cam_bg_B;

bg(2,1)=median(nu_intensity(untrans_id,2)./nu_intensity(untrans_id,1));
bg(2,2)=median(nu_intensity(untrans_id,3)./nu_intensity(untrans_id,1));
bg(2,3)=median(nu_intensity(untrans_id,4)./nu_intensity(untrans_id,1));

bg(3,1)=median((cell_intensity(untrans_id,2)-nu_intensity(untrans_id,2))./(cell_intensity(untrans_id,1)-nu_intensity(untrans_id,1)));
bg(3,2)=median((cell_intensity(untrans_id,3)-nu_intensity(untrans_id,3))./(cell_intensity(untrans_id,1)-nu_intensity(untrans_id,1)));
bg(3,3)=median((cell_intensity(untrans_id,4)-nu_intensity(untrans_id,4))./(cell_intensity(untrans_id,1)-nu_intensity(untrans_id,1)));

save([untransname filesep consname '_untrans.mat'],'untrans','bg')
save([transname filesep consname '_trans.mat'],'trans','bg')

%% Visusalization
% This section is designed to generate indicative images,
% allowing the selected transfected cells to be highlighted and visualized in the original images.
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
            % getIndex(z, c, t)  t-th timepoints in the z-th z-stacks in c-th
            % note getIndex starts from 0 but --- iPlane starts from 1
            I_R(:,:,ii) = bfGetPlane(reader, iPlane_R);
            I_G(:,:,ii) = bfGetPlane(reader, iPlane_G);
            I_B(:,:,ii) = bfGetPlane(reader, iPlane_B);
            imDAPI(:,:,ii) = bfGetPlane(reader, iPlane_DAPI);
        end

        I_raw=zeros(1024,1024,3,'uint16');
        I_raw(:,:,1) = max(I_R, [], 3);
        I_raw(:,:,3) = max(I_B, [], 3);

        I_raw_ad=I_raw;
        I_raw_ad(:,:,1)=imadjust(I_raw(:,:,1),[min(min(double(I_raw(:,:,1))))/65536,500/65536],[]);
        I_raw_ad(:,:,3)=autocontrast(I_raw(:,:,3),0.001,0.995);

        if isempty(trans(:,2)==num_of_image)
            continue
        end
        if isempty(untrans(:,2)==num_of_image)
            continue
        end
        trans_id_image=trans(trans(:,2)==num_of_image,1);
        untrans_id_image=untrans(untrans(:,2)==num_of_image,1);
        mask=imread(pathTomask);

        figure
        imshow(I_raw_ad)
        hold on
        for jjj=untrans_id_image'
            visboundaries((mask==jjj),'LineWidth',1,'LineStyle',':')
            hold on
        end
        saveas(gcf,[untransname filesep num2str(num_of_image) '_untrans.png'])

        figure
        imshow(I_raw_ad)
        hold on
        for jjj=trans_id_image'
            visboundaries((mask==jjj),'LineWidth',1,'LineStyle',':')
            hold on
        end
        saveas(gcf,[transname filesep num2str(num_of_image) '_trans.png'])
        close all
    end
end