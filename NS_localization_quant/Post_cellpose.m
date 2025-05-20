% Load excel sheets and add Bio-Formats library to MATLAB path
Fp=readcell(fullfile(pwd, 'current_constructs.xlsx'));
ExpLib=readtable(fullfile(pwd,'summary_cellpose_all.xlsx'));
Global_folder=readtable(fullfile(pwd,'global_folders.xlsx'));
addpath(fullfile(pwd,'bfmatlab'))
% Retrieve global folder paths
rawimage_folder=char(Global_folder.rawimage(1));
downstream_folder=char(Global_folder.downstream(1));
outline_folder=char(Global_folder.outline(1));
summary_folder=char(Global_folder.summary(1));
for m_file=1:size(Fp,1)
    consname = Fp{m_file,1};
    matchingRows = strcmp(consname, ExpLib.all_consname);
    [matched_id, ~] = find(matchingRows);
    if isempty(matched_id)
        errorMessage = 'Cannot match consname!';
        error(errorMessage);
    end
    folder = char(ExpLib.all_foldrer(matched_id));
    subfolder = char(ExpLib.all_subfolder(matched_id));
    % image parameters
    image_type = char(ExpLib.all_image_type(matched_id));
    image_name = char(ExpLib.all_image_name(matched_id));
    R_channel = ExpLib.all_R(matched_id);
    G_channel = ExpLib.all_G(matched_id);
    B_channel = ExpLib.all_B(matched_id);
    DAPI_channel = ExpLib.all_DAPI(matched_id);
    % Global image-processing quality control parameters
    nucleolar_removal = ExpLib.all_nucleolar_removal(matched_id);
    nucleolus_minV = ExpLib.all_nucleolus_minV(matched_id); 
    boundary_removal_cyto = ExpLib.all_boundary_removal_cyto(matched_id);
    boundary_removal_nu = ExpLib.all_boundary_removal_nu(matched_id);
    nucleus_minV = ExpLib.all_nucleus_minV(matched_id); %in pixel
    sat_num_threshold_R = ExpLib.all_sat_num_threshold_R(matched_id);
    sat_num_threshold_G = ExpLib.all_sat_num_threshold_G(matched_id);
    % threshold to select tranfected cells
    trans_threshold = ExpLib.all_trans_threshold(matched_id);
    % parameters for speckle segmentation
    num_of_std_B = ExpLib.all_num_of_std_B(matched_id);
    sigma1_B = ExpLib.all_sigma1_B(matched_id);
    sigma2_B = ExpLib.all_sigma2_B(matched_id);
    minvol_B =ExpLib.all_minvol_B(matched_id);
    minnum_B = ExpLib.all_minnum_B(matched_id);

    % create and define folders
    if isnumeric(folder)
        folder=num2str(folder);
    end
    if isnumeric(subfolder)
        folder=num2str(subfolder);
    end
    selectedimage_folder=[downstream_folder filesep folder filesep subfolder];
    create_folder
    cellname=[selectedimage_folder filesep 'cell'];
    nucleusname=[selectedimage_folder filesep 'nucleus'];
    untransname=[selectedimage_folder filesep 'untrans'];
    transname=[selectedimage_folder filesep 'trans'];
    spec_untrans_name=[selectedimage_folder filesep 'untrans' filesep 'spec'];
    spec_trans_name=[selectedimage_folder filesep 'trans' filesep 'spec'];
    %rna_trans_name=[selectedimage_folder filesep 'trans' filesep 'RNA'];

    disp(['post_cellpose:' num2str(m_file) '  ,  consname:' consname ';'])
    LoadCellBoundaries
    Auto_nucleus_Otsu
    SortTransfectedCells
    SegSpec_untrans
    SegSpec_trans

    remove_bg_R=ExpLib.all_remove_bg_R(matched_id);
    remove_bg_G=ExpLib.all_remove_bg_G(matched_id);
    remove_bg_B=ExpLib.all_remove_bg_B(matched_id);
    Spec_P


end

