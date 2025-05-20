%% Cellpose_input.m
%  This code is designed to generate input images suitable for cellpose
%  pre-trained models. Bio-Formats is used for reading raw images in nd2.

%  create folders to save the input images for cellpose
save_folder=[outline_folder filesep folder filesep subfolder];
if not(isfolder([save_folder]))
    mkdir([save_folder])
end

%%
if strcmp(image_type,'auto')
    pathToFile = [rawimage_folder filesep folder filesep subfolder filesep image_name '.nd2'];

    % From Bio-Formats importer:
    % the speedup gained from memoization will only happen after the first initialization of the reader for a particular file.
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader);
    reader.setId(pathToFile);

    % From Bio-Formats importer:
    % - getIndex(z, c, t): Retrieves the index of the t-th timepoint in the z-th Z-stack of the c-th channel.
    %   Note: Indexing starts from 0, not 1.
    % - int getSeriesCount(): Returns the total number of series in the file.

    seriesCount = reader.getSeriesCount();
    zstackCount=reader.getSizeZ();
    for i = 1:seriesCount
        DAPI_all_zstacks=uint16(zeros(1024,1024,zstackCount));
        I_R=uint16(zeros(1024,1024,zstackCount));
        reader.setSeries(i - 1); % Note: setSeries starts from 0
        for ii=1:zstackCount
            iPlane_R = reader.getIndex(ii-1, R_channel-1, 1-1) + 1;
            iPlane_DAPI = reader.getIndex(ii-1, DAPI_channel-1, 1-1) + 1;
            % Note: getIndex starts from 0 but --- iPlane starts from 1
            I_R(:,:,ii) = bfGetPlane(reader, iPlane_R);
            DAPI_all_zstacks(:,:,ii) = bfGetPlane(reader, iPlane_DAPI);
        end

        % This section of the code performs a maximum projection on the DAPI images.
        % A maximum projection combines the Z-stack images by taking the maximum intensity
        % value for each pixel across all Z-slices, effectively creating a 2D representation
        % that highlights the brightest regions.
        % After the projection, an auto-contrast adjustment is applied to enhance the image.
        DAPI_maxproj = max(DAPI_all_zstacks, [], 3);
        IR_maxproj = max(I_R, [], 3);
        Imin=min(min(DAPI_maxproj));
        Imax=max(max(DAPI_maxproj));
        lower=double(Imin);
        higher=lower+(double(Imax)-double(Imin))*0.3;
        imDAPI_ad = imadjust(DAPI_maxproj,[lower/65536 higher/65536],[],0.7);
        Imask = imDAPI_ad;

        mask_name=[save_folder filesep 'ad_' num2str(i) '.tif'];
        imwrite(Imask,mask_name)
    end

    reader.close()
end

%%
if strcmp(image_type,'manual')
    FileList = dir([rawimage_folder filesep folder filesep subfolder]);
    num_of_image=0;

    for jj=1:length(FileList)
        if ~FileList(jj).isdir
            num_of_image=num_of_image+1;
            pathToFile=[rawimage_folder filesep folder filesep subfolder filesep FileList(jj).name];
            reader = bfGetReader(pathToFile);

            seriesCount = reader.getSeriesCount();
            zstackCount=reader.getSizeZ();
            for i = 1:seriesCount
                DAPI_all_zstacks=uint16(zeros(1024,1024,zstackCount));
                I_R=uint16(zeros(1024,1024,zstackCount));
                reader.setSeries(i - 1); % Note: setSeries starts from 0
                for ii=1:zstackCount
                    iPlane_R = reader.getIndex(ii-1, R_channel-1, 1-1) + 1;
                    iPlane_DAPI = reader.getIndex(ii-1, DAPI_channel-1, 1-1) + 1;
                    % Note: getIndex starts from 0 but --- iPlane starts from 1
                    I_R(:,:,ii) = bfGetPlane(reader, iPlane_R);
                    DAPI_all_zstacks(:,:,ii) = bfGetPlane(reader, iPlane_DAPI);
                end
                DAPI_maxproj = max(DAPI_all_zstacks, [], 3);
                IR_maxproj = max(I_R, [], 3);
                Imin=min(min(DAPI_maxproj));
                Imax=max(max(DAPI_maxproj));
                lower=double(Imin);
                higher=lower+(double(Imax)-double(Imin))*0.3;
                imDAPI_ad = imadjust(DAPI_maxproj,[lower/65536 higher/65536],[],0.7);
                Imask = imDAPI_ad;


                mask_name=[save_folder filesep 'ad_' num2str(num_of_image) '.tif'];
                imwrite(Imask,mask_name)
            end

            reader.close()
        end
    end
end