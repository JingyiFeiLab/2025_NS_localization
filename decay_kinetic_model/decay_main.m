%% Global parameters
Fp=readcell('well-info.xlsx');
Global_folder=readtable(fullfile(pwd,'global_folder.xlsx'));
addpath(fullfile(pwd,'bfmatlab'))
rawimage_folder=char(Global_folder.rawimage(1));
outline_folder=char(Global_folder.outline(1));
adimage_folder=char(Global_folder.adimage(1));
selected_folder=char(Global_folder.selected(1));
if not(isfolder([selected_folder]))
    mkdir([selected_folder])
end
save_folder=char(Global_folder.save(1));
if not(isfolder([save_folder]))
    mkdir([save_folder])
end
summary_folder=char(Global_folder.summary(1));
if not(isfolder([summary_folder]))
    mkdir([summary_folder])
end
% Probe calibration is not required for fitting decay rates.
nor=1;
save([save_folder filesep 'nor_1.mat'],'nor')
save([save_folder filesep 'nor_2.mat'],'nor')
%% Extract single-cell intensities
NumberofCell = cell(size(Fp,1),9);
%NumberofCell:columns:
%1:wellname 2:# of raw cells 3:# of transfected cells 
% 4: median(IR) 5:median(IG) 
% 6 and 7: store reference thresholds (MinCRc, MinCGc) based on the
% robust median absolute deviation (MAD) method.
% 8 and 9: store the actual thresholds used for choosing transfected cells.
for m_file=1:size(Fp,1)
    disp(['Extraction:' num2str(m_file)])
    wellname=Fp{m_file,1};
    wellcons=Fp{m_file,2};
    probe_set=Fp{m_file,3};
    LoadCellBoundaries
    AutoNucleusOtsu
    NumberofCell{m_file,1}=wellcons;
    load([save_folder filesep wellcons 'cell_raw_ins.mat'])
    NumberofCell{m_file,2}=size(cell_raw_ins,1);
    cR = cell_raw_ins(:,12);
    cG = cell_raw_ins(:,13);
    NumberofCell{m_file,4}=median(cR);
    NumberofCell{m_file,5}=median(cG);
    MinCRc=5*1.482*mad(cR,1);
    MinCGc=5*1.482*mad(cG,1);
    NumberofCell{m_file,6}=MinCRc;
    NumberofCell{m_file,7}=MinCGc;
    NumberofCell{m_file,8}=median(cR)*1.15;
    NumberofCell{m_file,9}=median(cG)*1.15; % actual threshold_G
end

%% Background subtraction and selection of transfected cells
for m_file=1:size(Fp,1)
    disp(['Selection:' num2str(m_file)])
    wellname=Fp{m_file,1};
    wellcons=Fp{m_file,2};
    probe_set=Fp{m_file,3};
    SelectTransfected
    SubtractBackground
    % if m_file==1 || m_file==4 || m_file== 96 || m_file== 93
    %     QualityCheck
    % end
    if m_file==5 || m_file==8 || m_file== 1 || m_file== 4
        QualityCheck
    end
end
% QualityCheck ：This script is used to visualize the cells that were classified as 
% transfected cells. Specifically, it displays the original images with overlays or masks to 
% highlight selected transfected cells, enabling manual inspection of 
% selection accuracy.
%% Generate final statistics
for m_file=1:size(Fp,1)
    disp(['Statistics:' num2str(m_file)])
    wellname=Fp{m_file,1};
    wellcons=Fp{m_file,2};
    probe_set=Fp{m_file,3};
    load([save_folder filesep wellcons 'cell_raw_ins.mat'])
    NumberofCell{m_file,3}=size(cell_nobg_filtered,1);
end
save([save_folder filesep 'NumberofCell.mat'],'NumberofCell')
SummaryAllwells
SummaryAllcons_mean
% SummaryAllwells: summarizes the pre-mRNA intensities across all wells
% at each time point.
% SummaryAllcons_mean : performs fitting to determine 
% the RNA decay rate from time-course data.