Fp=readcell('well-info.xlsx');
Global_folder=readtable(fullfile(pwd,'global_folder.xlsx'));
addpath(fullfile(pwd,'bfmatlab'))
rawimage_folder=char(Global_folder.rawimage(1));
outline_folder=char(Global_folder.outline(1));
if not(isfolder([outline_folder]))
    mkdir([outline_folder])
end
%global parameters
for m_file=1:size(Fp,1)
    disp(['Loading well:' num2str(m_file)])
    if isnumeric(Fp{m_file,1})
        Fp{m_file,1}=num2str(Fp{m_file,1});
    end

    pathToFile = [rawimage_folder filesep Fp{m_file,1} '.nd2'];

    % the speedup gained from memoization will only happen after the first initialization of the reader for a particular file.
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader);
    reader.setId(pathToFile);

    % getIndex(z, c, t)  t-th timepoints in the z-th z-stacks in c-th
    % channels (note!!! starting from 0)
    % int getSeriesCount() Gets the number of series in this file.
    % reader = bfGetReader(pathToFile);
    seriesCount = reader.getSeriesCount();
    for i = 1:seriesCount
        reader.setSeries(i - 1); % note setSeries starts from 0
        %iPlane = reader.getIndex(iZ-1, iC -1, iT-1) + 1;
        % note getIndex starts from 0 but --- iPlane starts from 1
        I_raw = uint16(zeros(1024,1024,3));
        I_raw(:,:,1) = bfGetPlane(reader, 1);
        I_raw(:,:,2) = bfGetPlane(reader, 2);
        I_raw(:,:,3) = bfGetPlane(reader, 3);
        [m,n]=size(I_raw(:,:,1));
        I_ad = I_raw;
        % % Channel assignment for imaging:
        %
        % Channel 1: DAPI 
        % Channel 2: CF568 (for exonic probes)
        % Channel 3: AF647 (for intronic probes)
        I_ad(:,:,3) = autocontrast(I_raw(:,:,1),(50)/(m*n),1-(50)/(m*n));
        I_ad(:,:,2) = autocontrast(I_raw(:,:,2),(50)/(m*n),1-(50)/(m*n));
        I_ad(:,:,1) = autocontrast(I_raw(:,:,3),(50)/(m*n),1-(50)/(m*n));
        acdc_name=[outline_folder filesep 'ad-' Fp{m_file,1} '_' num2str(i) '.tif'];
        imwrite(I_ad,acdc_name)
    end

    reader.close()


    %     Iseries=bfopen([rawimage_folder filesep Fp{m_file,1} '.nd2']);
    %
    %     for i=1:size(Iseries,1)
    %         I_raw = uint16(zeros(1024,1024,3));
    %         I_raw(:,:,1) = Iseries{i,1}{1,1};
    %         I_raw(:,:,2) = Iseries{i,1}{2,1};
    %         I_raw(:,:,3) = Iseries{i,1}{3,1};
    %         [m,n]=size(I_raw(:,:,1));
    %         I_ad = I_raw;
    %         I_ad(:,:,1) = autocontrast(I_raw(:,:,1),(50)/(m*n),1-(50)/(m*n));
    %         I_ad(:,:,2) = autocontrast(I_raw(:,:,2),(50)/(m*n),1-(50)/(m*n));
    %         I_ad(:,:,3) = autocontrast(I_raw(:,:,3),(50)/(m*n),1-(50)/(m*n));
    %         acdc_name=[save_folder filesep 'ad-' Fp{m_file,1} '_' num2str(i) '.tif'];
    %         imwrite(I_ad,acdc_name)
    %     end
end
