% Load data from previous steps
load([nucleusname filesep consname 'nu_raw.mat'])
load([transname filesep consname '_trans.mat'])
trans_id=trans(:,3);
trans_num=length(trans_id);
% The matrix 'spec_intensity' consists of multiple columns,
% where each column represents a specific characteristic (at the single-cell level):
% - col1: 'Area' - Represents the total area of the speckles.
% - col2:  - Indicates the intensity value inside speckles for the first channel.
% - col3:  - Represents the intensity value inside speckles for the second channel.
% - col4:  - Represents the intensity value for the third channel.
% - col5:  - Contains cell ID (the identifier in nu_raw.mat).
% - col6:  - Specifies the number of speckles per cell.
spec_intensity=zeros(trans_num,6);
% spec_prop: A structure containing quantitative properties of detected
% speckles (at the single-speckle level).
spec_prop = struct();
% j is used to count cells that have more than minnum_B speckles.
j=0;
for i = 1:trans_num
    id=trans_id(i);
    % Extract I_crop and bw_dapi from "nu_raw.mat"
    % Specify channel to be Alexa Fluor 488
    I_crop=nu_raw{id,1};
    if nucleolar_removal == 1
        bw_dapi=nu_raw{id,3};
    else
        bw_dapi=nu_raw{id,2};
    end
    channel=3;
    % speckle segmentation
    % quality check: cells with speckles less than minnum_B are removed
    [speckleprop,bw_spec] = speckle_thre(bw_dapi,I_crop,channel,minvol_B,sigma1_B,sigma2_B,num_of_std_B,id);
    if speckleprop.N < minnum_B
        speckleprop=[];
        bw_spec=[];
        save([spec_trans_name filesep 'speckle_', num2str(i), '.mat'],'speckleprop','bw_spec');
        continue
    end

    j=j+1;
    % To visualize speckle segmentation, I_crop is enhanced and resized
    if j<10
        scale=4;
        I_raw_ad=I_crop;
        I_crop_ad=I_crop.*uint16(bw_dapi);
        [m,n]=size(I_raw_ad);
        I_raw_ad(:,:,1)=autocontrast(I_crop_ad(:,:,1)+uint16(~bw_dapi).*min(min(nonzeros(I_crop_ad(:,:,1)))),5/(m*n),1-5/(m*n));
        I_raw_ad(:,:,2)=autocontrast(I_crop_ad(:,:,2)+uint16(~bw_dapi).*min(min(nonzeros(I_crop_ad(:,:,2)))),5/(m*n),1-5/(m*n));
        I_raw_ad(:,:,3)=autocontrast(I_crop_ad(:,:,3)+uint16(~bw_dapi).*min(min(nonzeros(I_crop_ad(:,:,3)))),5/(m*n),1-5/(m*n));
        I_raw_mag=imresize(I_raw_ad,scale);
        bw_spec_mag=imresize(bw_spec,scale);
        figure
        imshow(I_raw_mag(:,:,3))
        hold on
        visboundaries(bw_spec_mag,'LineWidth',0.5,'LineStyle',':')
        saveas(gcf,[spec_trans_name filesep num2str(j) '_specB.png'])
        figure
        imshow(I_raw_mag)
        hold on
        visboundaries(bw_spec_mag,'LineWidth',0.5,'LineStyle',':')
        saveas(gcf,[spec_trans_name filesep num2str(j) '_specRGB.png'])
        imwrite(I_crop,[spec_trans_name filesep num2str(j) '.tif'])
        imwrite(bw_spec,[spec_trans_name filesep num2str(j) '_spec_bw.png'])

    end

    save([spec_trans_name filesep 'speckle_', num2str(i), '.mat'],'speckleprop','bw_spec');
    spec_prop(j).singlecell = speckleprop;
    spec_intensity(j,2)=sum(sum(uint16(bw_spec).*I_crop(:,:,1)));
    spec_intensity(j,3)=sum(sum(uint16(bw_spec).*I_crop(:,:,2)));
    spec_intensity(j,4)=sum(sum(uint16(bw_spec).*I_crop(:,:,3)));
    spec_intensity(j,1)=sum(sum(uint16(bw_spec)));
    spec_intensity(j,5)=id;
    spec_intensity(j,6)=speckleprop.N;
end
close all
spec_intensity(spec_intensity(:,1)==0,:)=[];
save([spec_trans_name filesep consname 'spec_intensity.mat'],'spec_intensity','spec_prop')

