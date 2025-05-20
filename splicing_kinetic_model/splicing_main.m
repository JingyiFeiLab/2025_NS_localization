%% Global path and well-level info
% This section initializes file paths and folders required for downstream 
% image analysis and data processing.
Fp=readcell('well-info.xlsx');
Global_folder=readtable(fullfile(pwd,'global_folder.xlsx'));
addpath(fullfile(pwd,'bfmatlab'))
rawimage_folder=char(Global_folder.rawimage(1));
outline_folder=char(Global_folder.outline(1));
allcells_folder=char(Global_folder.allcells(1));
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
% The first or second well: (Plad B treated)
% to be used for calibrating the probes
m_file_cal_1=1; % or enter the well for calibration
m_file_cal_2=0;
%% Single-cell statistics
% This section prepares a summary matrix `NumberofCell` to extract and store 
% cell-level statistics for each well defined in `Fp`. 
%NumberofCell:columns
% 1:wellname 2:# of raw cells 
% 3:# of transfected cells 
% 4 and 5: Background estimated by non-cell regions
% 8 and 9: store reference thresholds (MinCRc, MinCGc) based on the
% robust median absolute deviation (MAD) method.
% 6 and 7: store the actual thresholds used for choosing transfected cells.
for m_file=1:size(Fp)
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
% col 4 and 5 : Background estimated by non-cell regions
EstimateBgFromAllCells
load([save_folder filesep wellcons 'bg_from_allcells.mat'])
%NumberofCell{m_file,4}=bg_info(:,1);
%NumberofCell{m_file,5}=bg_info(:,2);
NumberofCell{m_file,4}=bg_R_fromallcells;
NumberofCell{m_file,5}=bg_G_fromallcells;
MinCRc=5*1.482*mad(cR,1);
MinCGc=5*1.482*mad(cG,1);
% `MinCRc` and `MinCGc` gives a reference value of intensity threshold chose for 
% distinguishing transfected cells.
% It approximates 5 standard deviations above the median under the assumption 
% of a normal distribution, where 1.482 is a scaling factor converting MAD 
% to standard deviation.
NumberofCell{m_file,6}=NumberofCell{m_file,4}*5/130; % actual threshold_R
NumberofCell{m_file,7}=NumberofCell{m_file,5}*5/130; % actual threshold_G
NumberofCell{m_file,8}=MinCRc;
NumberofCell{m_file,9}=MinCGc;
end
%% Background Correction and Transfected Cell Selection
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
% QualityCheck ：This script is used to visualize the cells that were classified as 
% transfected cells. Specifically, it displays the original images with overlays or masks to 
% highlight selected transfected cells, enabling manual inspection of 
% selection accuracy.
if m_file==1 || m_file==6 || m_file== 7 || m_file== 13
     QualityCheck
end
end
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
% Probe Calibration
disp('ProbeCalibration')
wellname=Fp{m_file_cal_1,1};
wellcons=Fp{m_file_cal_1,2};
probe_set=Fp{m_file_cal_1,3};
% This script estimates the relative fluorescence efficiency between 
% intronic and exonic probes by comparing their signal intensities in 
% a calibration dataset. The goal is to determine the correction factor 
% needed to normalize their signals for downstream kinetic modeling.
% In the imaging setup used, probe concentrations and imaging parameters 
% were pre-adjusted so that the efficiency ratio between intronic and 
% exonic probes is approximately 1. Therefore, the fitted calibration 
% result is expected to yield a ratio close to 1.
ProbesCalibration
% wellname=Fp{m_file_cal_2,1};
% wellcons=Fp{m_file_cal_2,2};
% probe_set=Fp{m_file_cal_2,3};
% ProbesCalibration
add_info=readcell('additional_info.xlsx');
add_info=add_info(2:end,:); % remove the well used for probe calibration
consname_info=add_info(:,1);
mu=cell2mat(add_info(:,3));
delta_p=cell2mat(add_info(:,4));
% SummaryAllwells_splicing : summarizes the splicing data across all wells by extracting 
% the relative abundance of pre-mRNA ([P]) and spliced RNA ([S]) at each 
% time point.
SummaryAllwells_splicing
% SummaryAllcons_splicing : performs kinetic fitting of splicing data for each individual construct.
SummaryAllcons_splicing