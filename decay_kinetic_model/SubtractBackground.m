load([save_folder filesep wellcons 'cell_raw_ins.mat'])
cell_num_all=size(cell_filtered,1);
cell_nobg_filtered=cell_filtered;
% bg_R=median(cell2mat(cell_raw(:,15)));
% bg_G=median(cell2mat(cell_raw(:,16)));
bg_R = NumberofCell{m_file,4};
bg_G = NumberofCell{m_file,5};
bg_info(:,1)=bg_R;
bg_info(:,2)=bg_G;
for i=1:cell_num_all
    A_cell=cell_filtered(i,10);
    A_nu=cell_filtered(i,11);
    cell_nobg_filtered(i,4)=cell_filtered(i,4)-A_cell*bg_R;
    cell_nobg_filtered(i,5)=cell_filtered(i,5)-A_cell*bg_G;
    cell_nobg_filtered(i,6)=cell_filtered(i,6)-A_nu*bg_R;
    cell_nobg_filtered(i,7)=cell_filtered(i,7)-A_nu*bg_G;
    cell_nobg_filtered(i,8)=cell_filtered(i,8)-(A_cell-A_nu).*bg_R;
    cell_nobg_filtered(i,9)=cell_filtered(i,9)-(A_cell-A_nu).*bg_G;

    cell_nobg_filtered(i,12)=cell_filtered(i,12)-bg_R;
    cell_nobg_filtered(i,13)=cell_filtered(i,13)-bg_G;
    cell_nobg_filtered(i,14)=cell_filtered(i,14)-bg_R;
    cell_nobg_filtered(i,15)=cell_filtered(i,15)-bg_G;
    if cell_filtered(i,16)==0 % if there is no cytoplasmic signal, don't subtract background
    cell_nobg_filtered(i,16)=cell_filtered(i,16);
    else
    cell_nobg_filtered(i,16)=cell_filtered(i,16)-bg_R;
    end
    if cell_filtered(i,17)==0 % if there is no cytoplasmic signal, don't subtract background
    cell_nobg_filtered(i,17)=cell_filterd(i,17);
    else
    cell_nobg_filtered(i,17)=cell_filtered(i,17)-bg_G;
    end

end
save([save_folder filesep wellcons 'cell_raw_ins.mat'],'cell_raw_ins','bg_info','cell_nobg_filtered','cell_filtered')