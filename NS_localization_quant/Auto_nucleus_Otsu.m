%% Using Otsu's method to segment nucleus
% Outputs descriptions:
% 1)'nu_intensity' is a matrix with the following columns:
% - col1: Area: The area of the nucleus.
% - col2: total nuclear intensity of the first channel.
% - col3: total nuclear intensity of the second channel.
% - col4: total nuclear intensity of the third channel.
% 2) 'nu_raw' is a cell array with the following columns:
% - col1: 'I_raw_crop' - The cropped raw image data for the cell.
% - col2: 'bw_nucleus' - Binary mask of the segmented nucleus.
% - col3: 'bw_nucleus_r' - Binary mask of the segmented nucleus with
% nucleoli removed.
% - col4: 'imDAPI_crop' - Cropped DAPI image.

%%
load([cellname filesep consname 'cell_raw.mat'])
% Pre-allocate memory for variables
cell_num_all=size(cell_raw,1);
nu_intensity=zeros(cell_num_all,4);
nu_intensity_r=zeros(cell_num_all,4);
nu_raw=cell(cell_num_all,4);

for i=1:cell_num_all
    I_raw=cell_raw{i,1};
    bw_cyto=cell_raw{i,2};
    [m,n]=size(I_raw(:,:,1));

    imDAPI=autocontrast(cell_raw{i,3},0.001,0.999);
    imDAPI_ad=imDAPI.*uint16(bw_cyto)+min(nonzeros(imDAPI.*uint16(bw_cyto))).*(uint16(~bw_cyto));

    % Segment nucleus by Otsu's method
    I_thre=imbinarize(imDAPI_ad);
    seD = strel('diamond',1);
    Ibw = imdilate(I_thre,seD);
    Ibw = imfill(Ibw,'holes');
    Ibw = imerode(Ibw,seD);
    Ibw = bwareaopen(Ibw,nucleus_minV);
    bw_nucleus= bw_cyto & Ibw;

    % Quality check : if no nucleus detected, continue
    if isempty(nonzeros(uint16(bw_nucleus)))
        continue
    end
    % Boundary check :
    % boundary_removal_cyto==1 : If the cell's edge touches the image boundary, the cell is removed.
    % boundary_removal_nu==1 : If the cell's nucleus touches the image boundary, the cell is removed.
    if boundary_removal_cyto==1 && cell_raw{i,8}==1
        continue
    end
    [rows, cols] = find(bw_nucleus);
    [height, width] = size(bw_nucleus);
    touchesEdge_nu = any(rows == 1) || any(rows == height) || any(cols == 1) || any(cols == width);
    if (boundary_removal_nu==1) && touchesEdge_nu
        continue
    end

    % This part is to remove nucleoli(holes in DAPI nucleus images)
    if nucleolar_removal == 1
        I_comb=I_raw(:,:,3)+cell_raw{i,3};
        I_pseudo=min(nonzeros(uint16(bw_nucleus).*I_comb)).*(uint16(ones(size(I_comb)))-uint16(bw_nucleus))+I_comb.*uint16(bw_nucleus);
        I_pseudo_ad=autocontrast(I_pseudo,0.05,0.95);

        bw_nucleus_r=imbinarize(I_pseudo_ad);
        bw_nucleolus=(~bw_nucleus_r);
        bw_nucleolus = bwareaopen(bw_nucleolus,nucleolus_minV);
        bw_nucleus_r=(~bw_nucleolus);
        bw_nucleus_filled=imfill(bw_nucleus_r,'holes');
        bw_nucleus_boundary=(~bw_nucleus_filled)&(bw_nucleus);
        bw_nucleus_r=(bw_nucleus_r)|(bw_nucleus_boundary);
    else
        bw_nucleus_r=bw_nucleus;
    end

    if i<10
        figure
        I_ad=I_raw;
        I_ad(:,:,1)=imDAPI_ad;
        I_ad(:,:,2)=imDAPI_ad;
        I_ad(:,:,3)=autocontrast(I_raw(:,:,3),0.001,0.995);
        imshow(I_ad)
        hold on
        visboundaries(bw_nucleus_r,'LineWidth',0.5,'LineStyle',':')
        saveas(gcf,[nucleusname filesep num2str(i) '_nucleus_removed.png'])
        figure
        imshow(I_ad)
        hold on
        visboundaries(bw_nucleus,'LineWidth',0.5,'LineStyle',':')
        saveas(gcf,[nucleusname filesep num2str(i) '_nucleus.png'])
    end


    nu_intensity(i,2)=sum(sum(uint16(bw_nucleus).*I_raw(:,:,1)));
    nu_intensity(i,3)=sum(sum(uint16(bw_nucleus).*I_raw(:,:,2)));
    nu_intensity(i,4)=sum(sum(uint16(bw_nucleus).*I_raw(:,:,3)));
    nu_intensity(i,1)=sum(sum(uint16(bw_nucleus)));

    nu_intensity_r(i,2)=sum(sum(uint16(bw_nucleus_r).*I_raw(:,:,1)));
    nu_intensity_r(i,3)=sum(sum(uint16(bw_nucleus_r).*I_raw(:,:,2)));
    nu_intensity_r(i,4)=sum(sum(uint16(bw_nucleus_r).*I_raw(:,:,3)));
    nu_intensity_r(i,1)=sum(sum(uint16(bw_nucleus_r)));

    nu_raw{i,1}=I_raw;
    nu_raw{i,2}=bw_nucleus;
    nu_raw{i,3}=bw_nucleus_r;
    nu_raw{i,4}=cell_raw{i,3};
end
save([nucleusname filesep consname 'nu_intensity.mat'],'nu_intensity','nu_intensity_r')
save([nucleusname filesep consname 'nu_raw.mat'],'nu_raw')
close all