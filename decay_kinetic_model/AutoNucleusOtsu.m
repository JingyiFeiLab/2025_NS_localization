load([save_folder filesep wellcons 'cell_raw.mat'])
load([save_folder filesep wellcons 'cell_raw_ins.mat'])
cell_num_all=size(cell_raw,1);
for i=1:cell_num_all
    I_raw=cell_raw{i,1};
    bw_cyto=cell_raw{i,2};
    [m,n]=size(I_raw(:,:,1));

    imDAPI=I_raw(:,:,3);
    imDAPI_mer=imDAPI.*uint16(bw_cyto)+min(nonzeros(imDAPI.*uint16(bw_cyto))).*(uint16(~bw_cyto));
    imDAPI_ad=autocontrast(imDAPI_mer,(5)/(m*n),1-(5)/(m*n));

    I_thre=imbinarize(imDAPI_ad);
    seD = strel('diamond',1);
    Ibw = imdilate(I_thre,seD);
    Ibw = imfill(Ibw,'holes');
    Ibw = imerode(Ibw,seD);
    %Ibw = bwareaopen(Ibw,nucleus_minV);
    bw_nucleus= bw_cyto & Ibw;
    %     if isempty(nonzeros(uint16(bw_nucleus)))
    %         continue
    %     end

    if i<2
        figure
        imshow(imDAPI_ad)
        hold on
        visboundaries(bw_nucleus,'LineWidth',0.5,'LineStyle',':')
        saveas(gcf,[save_folder filesep wellcons '_' num2str(i) '_nucleus.png'])
    end

    cell_raw{i,3}=bw_nucleus;
    cell_raw_ins(i,4)=sum(sum(uint16(bw_cyto).*I_raw(:,:,1)));
    cell_raw_ins(i,5)=sum(sum(uint16(bw_cyto).*I_raw(:,:,2)));
    cell_raw_ins(i,6)=sum(sum(uint16(bw_nucleus).*I_raw(:,:,1)));
    cell_raw_ins(i,7)=sum(sum(uint16(bw_nucleus).*I_raw(:,:,2)));
    cell_raw_ins(i,8)=cell_raw_ins(i,4)-cell_raw_ins(i,6);
    cell_raw_ins(i,9)=cell_raw_ins(i,5)-cell_raw_ins(i,7);
    A_cell=sum(sum(uint16(bw_cyto)));
    cell_raw_ins(i,10)=A_cell;
    A_nu=sum(sum(uint16(bw_nucleus)));
    cell_raw_ins(i,11)=A_nu;
    cell_raw_ins(i,12)=cell_raw_ins(i,4)./A_cell;
    cell_raw_ins(i,13)=cell_raw_ins(i,5)./A_cell;
    cell_raw_ins(i,14)=cell_raw_ins(i,6)./A_nu;
    cell_raw_ins(i,15)=cell_raw_ins(i,7)./A_nu;
    if A_cell == A_nu  % if A_cell=A_nu, this means cytoplasm is not properly indentified, this could be due to no cytoplamsic signals
        cell_raw_ins(i,16)=0;
        cell_raw_ins(i,17)=0;
    else
        cell_raw_ins(i,16)=cell_raw_ins(i,8)./(A_cell-A_nu);
        cell_raw_ins(i,17)=cell_raw_ins(i,9)./(A_cell-A_nu);
    end
end
save([save_folder filesep wellcons 'cell_raw.mat'],'cell_raw')
save([save_folder filesep wellcons 'cell_raw_ins.mat'],'cell_raw_ins','bg_info')
close all