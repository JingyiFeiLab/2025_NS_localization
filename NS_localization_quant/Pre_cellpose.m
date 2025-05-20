%% Pre_cellpose.m
%  This code is designed to generate input images suitable for cellpose
%  pre-trained models. It reads construct information, global folder paths, 
%  and experimental parameters from input Excel files.
%  Dependency:
%  excel files: current_constructs.xlsx, summary_cellpose_all.xlsx, and
%  global_folders.xlsx;
%  scripts: Cellpose_input.m
%% 

% Load input files
Fp = readcell(fullfile(pwd, 'current_constructs.xlsx')); % Load the list of constructs
ExpLib = readtable(fullfile(pwd, 'summary_cellpose_all.xlsx')); % Load experiment data summary
Global_folder = readtable(fullfile(pwd, 'global_folders.xlsx')); % Load global folder paths
addpath(fullfile(pwd, 'bfmatlab')) % Add Bio-Formats library to MATLAB path

% Retrieve global folder paths
rawimage_folder = char(Global_folder.rawimage(1)); % Path to raw images
downstream_folder = char(Global_folder.downstream(1)); % Path to downstream folder
outline_folder = char(Global_folder.outline(1)); % Path to outline folder
summary_folder = char(Global_folder.summary(1)); % Path to summary folder

%% Loop through each construct
for m_file = 1:size(Fp, 1)
    consname = Fp{m_file, 1}; % Get construct name
    matchingRows = strcmp(consname, ExpLib.all_consname); % Find matching rows in the experiment library
    [matched_id, ~] = find(matchingRows); % Get the indices of matched rows
    
    % Check if the construct is found in the experiment library
    if isempty(matched_id)
        errorMessage = 'Cannot match consname!'; % Error message for unmatched construct
        error(errorMessage); % Stop execution with an error
    end

    % Retrieve parameters for the matched construct
    folder = char(ExpLib.all_foldrer(matched_id)); % Get main folder name
    subfolder = char(ExpLib.all_subfolder(matched_id)); % Get subfolder name
    image_type = char(ExpLib.all_image_type(matched_id)); % Image type
    image_name = char(ExpLib.all_image_name(matched_id)); % Image name
    R_channel = ExpLib.all_R(matched_id); % R channel
    G_channel = ExpLib.all_G(matched_id); % G channel
    B_channel = ExpLib.all_B(matched_id); % B channel
    DAPI_channel = ExpLib.all_DAPI(matched_id); % DAPI channel

    % Handle numeric folder names by converting to strings if needed
    if isnumeric(folder)
        folder = num2str(folder);
    end
    if isnumeric(subfolder)
        folder = num2str(subfolder);
    end

    % Display progress for debugging or monitoring
    disp(['pre_cellpose:' num2str(m_file) '  ,  consname:' consname ';'])

    % Call the Cellpose input processing script or function
    Cellpose_input
end
