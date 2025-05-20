load([save_folder filesep wellcons 'cell_raw_ins.mat'])
%display selected cell info
Image_id=cell_nobg_filtered(:,2);
for i=1:max(Image_id)
    bw_selected_cell=false(1024);
    cell_in_image=cell_nobg_filtered(Image_id==i,3);
    pathTomask=[outline_folder filesep 'ad-' wellname '_' num2str(i) '_cp_masks.png'];
    if isempty(cell_in_image)
        continue
    end
    mask=imread(pathTomask);
    ad_image=imread([adimage_folder filesep 'ad-' wellname '_' num2str(i) '.tif']);
    ad_image(:,:,1)=imadjust(ad_image(:,:,1),[0 10000/65536],[]);
    ad_image(:,:,2)=imadjust(ad_image(:,:,2),[0 10000/65536],[]);
    for j=1:size(cell_in_image,1)
        cell_id_image=cell_in_image(j);
        bw_selected_cell = (bw_selected_cell | mask == cell_id_image);
    end
    figure
    imshow(ad_image)
    hold on
    visboundaries(bw_selected_cell,'LineWidth',0.5)
    ax = gca;
    exportgraphics(ax,[selected_folder filesep wellname '_' num2str(i) '_filtered.png'])
end
close all
% display original selected cells
Image_id=cell_raw_ins(:,2);
for i=1:max(Image_id)
    bw_selected_cell=false(1024);
    cell_in_image=cell_raw_ins(Image_id==i,3);
    pathTomask=[outline_folder filesep 'ad-' wellname '_' num2str(i) '_cp_masks.png'];
    if isempty(cell_in_image)
        continue
    end
    mask=imread(pathTomask);
    ad_image=imread([adimage_folder filesep 'ad-' wellname '_' num2str(i) '.tif']);
    ad_image(:,:,3)=imadjust(ad_image(:,:,3),[0 10000/65536],[]);
    for j=1:size(cell_in_image,1)
        cell_id_image=cell_in_image(j);
        bw_selected_cell = (bw_selected_cell | mask == cell_id_image);
    end
    figure
    imshow(ad_image)
    hold on
    visboundaries(bw_selected_cell,'LineWidth',0.5)
    ax = gca;
    exportgraphics(ax,[selected_folder filesep wellname '_' num2str(i) '_raw.png'])
end
close all

