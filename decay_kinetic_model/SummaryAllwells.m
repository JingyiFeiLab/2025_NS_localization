sum_allwell=cell(size(Fp,1),19);
%date consname time 
%4 5 6 mean(IR_total_cell(raw)) std se
%7 8 9 mean(IG_total_cell(raw)) std se
%10 11 12 IR_total_cell(filter) std se
%13 14 15 IG_total_cell(filter) std se
%16 17 #raw cells #filtered cells
%18 19 vector IR_total_filter IG_total_filter
load([save_folder filesep 'NumberofCell.mat'])
for m_file=1:size(Fp,1)
wellname=Fp{m_file,1};
wellcons=Fp{m_file,2};
probe_set=Fp{m_file,3};
if NumberofCell{m_file,3}==0
    continue;
end
load([save_folder filesep wellcons 'cell_raw_ins.mat'])

IRc_raw=cell_raw_ins(:,4);
IGc_raw=cell_raw_ins(:,5);
IRc_filtered=cell_nobg_filtered(:,4);
IGc_filtered=cell_nobg_filtered(:,5);

sum_allwell{m_file,4}=mean(IRc_raw);
sum_allwell{m_file,5}=std(IRc_raw);
sum_allwell{m_file,6}=std(IRc_raw)./sqrt(length(IRc_raw));
sum_allwell{m_file,7}=mean(IGc_raw);
sum_allwell{m_file,8}=std(IGc_raw);
sum_allwell{m_file,9}=std(IGc_raw)./sqrt(length(IGc_raw));
sum_allwell{m_file,10}=mean(IRc_filtered);
sum_allwell{m_file,11}=std(IRc_filtered);
sum_allwell{m_file,12}=std(IRc_filtered)./sqrt(length(IRc_filtered));
sum_allwell{m_file,13}=mean(IGc_filtered);
sum_allwell{m_file,14}=std(IGc_filtered);
sum_allwell{m_file,15}=std(IGc_filtered)./sqrt(length(IGc_filtered));
sum_allwell{m_file,18}=IRc_filtered;
sum_allwell{m_file,19}=IGc_filtered;

sum_allwell{m_file,16}=NumberofCell{m_file,2};
sum_allwell{m_file,17}=NumberofCell{m_file,3};
pattern = '(\d{8})-(\w+)-(\w+)-([\d.]+)h_';
tokens = regexp(wellcons, pattern, 'tokens');
sum_allwell{m_file,1} = tokens{1}{1};
sum_allwell{m_file,2} = [tokens{1}{2} '-' tokens{1}{3}];
sum_allwell{m_file,3} = tokens{1}{4};
end
sum_allwell=sum_allwell(~cellfun('isempty',sum_allwell(:,1)),:);
save([summary_folder filesep 'sum_allwell.mat'],'sum_allwell')