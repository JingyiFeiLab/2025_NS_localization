load([save_folder filesep wellcons 'cell_raw_ins.mat'])
IRn=cell_raw_ins(:,4);
IGn=cell_raw_ins(:,5);
ft = fittype({'x'});
[Cali,gof_Cali] = fit(IRn,IGn,ft);
nor = coeffvalues(Cali);
if probe_set == 1
    save([save_folder filesep 'nor_1.mat'],'nor','gof_Cali')
else
    save([save_folder filesep 'nor_2.mat'],'nor','gof_Cali')
end