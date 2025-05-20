load([save_folder filesep wellcons 'cell_raw_ins.mat'])
cell_filtered=cell_raw_ins;
cell_num_all=size(cell_raw_ins,1);

MinAnu=20;
% IR=cell2mat(cell_raw(:,15));
% IG=cell2mat(cell_raw(:,16));
% MinCRc=5*1.482*mad(IR,1);
% MinCGc=5*1.482*mad(IG,1);
bg_R = NumberofCell{m_file,4};
bg_G = NumberofCell{m_file,5};
bg_info(:,1)=bg_R;
bg_info(:,2)=bg_G;
% MinCRc = NumberofCell{m_file,6};
% MinCGc = NumberofCell{m_file,7};
thre_G=NumberofCell{m_file,9};
for i=1:cell_num_all

    % filter cells without nucleus
    if cell_filtered(i,11)<MinAnu
        cell_filtered(i,1)=0;
    end
    % filter cells without nuclear intronic signals
    if cell_filtered(i,13)<thre_G
        cell_filtered(i,1)=0;
    end

end
%cell_filtered=cell_filtered(~cellfun('isempty',cell_filtered(:,1)),:);
cell_filtered(cell_filtered(:,1)==0,:)=[];
save([save_folder filesep wellcons 'cell_raw_ins.mat'],'bg_info','cell_raw_ins','cell_filtered')