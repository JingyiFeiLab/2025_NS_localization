function [speckleprop,bw_spec] = speckle_thre(bw_dapi,I_crop,channel,minV,sigma1,sigma2,num_of_std,cell_id)
%% Input description
% bw_dapi(logical):a binary mask representing the segmented DAPI image
% I_crop(3-channel): the cropped image containing the region of interest.
% channel: to specify the channel to analyze
% minV: lower limit of speckle size
% sigma1, sigma2,num_of_std: used in DoG filter and threshold calcualtion
% cell_id: A unique identifier for the cell being processed.
%% Output description
% speckleprop(struct): A structure containing quantitative properties of detected speckles in the processed image.
% bw_spec(logical): A binary mask representing the detected speckles.
%% Normalization
% Normalization: rescale image channel for segmentation by a linear transformation
if channel == 1
    I_nu=uint16(bw_dapi).*I_crop(:,:,1);
end
if channel == 2
    I_nu=uint16(bw_dapi).*I_crop(:,:,2);
end
if channel == 3
    I_nu=uint16(bw_dapi).*I_crop(:,:,3);
end
I_mask_raw = I_nu+uint16(~bw_dapi).*min(min(nonzeros(I_nu)));
[m,n] = size(I_crop);
I_mask = autocontrast(I_mask_raw,5/(m*n),1-5/(m*n)).*uint16(bw_dapi);

I_R(:,:) = I_crop(:,:,1);
I_G(:,:) = I_crop(:,:,2);
I_B(:,:) = I_crop(:,:,3);

%% Thresholding
% First apply a DoG filter, then define the intensity threshold
I_sub=DoG(I_mask,sigma1,sigma2);
intensity_raw=reshape(double(bw_dapi).*double(I_sub),[],1);
intensity_thre=mean(intensity_raw)+num_of_std*std(intensity_raw);

Ibw = imbinarize(I_sub,intensity_thre/65536);  % thresholding
Ibw1 = imfill(Ibw,'holes');			% fills holes in the image
Ibw2 = bwareaopen(Ibw1,minV);       % remove tiny speckles
bw_spec = bw_dapi & Ibw2;


%% Find connected objects
% Find and count connected components in the obtained binary image
cc = bwconncomp(bw_spec, 8);
Pobj3D = regionprops(cc, 'Area','BoundingBox','Centroid','PixelList','PixelIdxList','Circularity','Eccentricity');
Pobj3D_R = regionprops(cc,I_R,"MeanIntensity");
Pobj3D_G = regionprops(cc,I_G,"MeanIntensity");
Pobj3D_B = regionprops(cc,I_B,"MeanIntensity");
bw_nucleoplasm = (bw_dapi) & (~bw_spec);
IR_nucleoplasm = sum(sum(uint16(bw_nucleoplasm).*I_R))./sum(sum(uint16(bw_nucleoplasm)));
IG_nucleoplasm = sum(sum(uint16(bw_nucleoplasm).*I_G))./sum(sum(uint16(bw_nucleoplasm)));
IB_nucleoplasm = sum(sum(uint16(bw_nucleoplasm).*I_B))./sum(sum(uint16(bw_nucleoplasm)));
N = length(Pobj3D);

%% Write Outputs
% speckleprop: cell_id, speckle_id
% Area, Circularity, Eccentricity,
% meanIR, meanIG, meanIB,
% IR_nucleoplasm, IG_nucleoplasm, IB_nucleoplasm
speckleprop.N = N;
speckleprop.cell_id = cell_id;
for i = 1:N
    % Note: In MATLAB's image coordinate system:
    % x represents the column index (horizontal position).
    % y represents the row index (vertical position).
    speckleprop.spec_id = i;
    speckleprop.col(i) = Pobj3D(i).Centroid(1);
    speckleprop.row(i) = Pobj3D(i).Centroid(2);
    speckleprop.v(i) = Pobj3D(i).Area;
    speckleprop.Circularity(i)=Pobj3D(i).Circularity;
    speckleprop.Eccentricity(i)=Pobj3D(i).Eccentricity;

    speckleprop.meanIR(i)=Pobj3D_R(i).MeanIntensity;
    speckleprop.meanIG(i)=Pobj3D_G(i).MeanIntensity;
    speckleprop.meanIB(i)=Pobj3D_B(i).MeanIntensity;
    speckleprop.IR_nucleoplasm=IR_nucleoplasm;
    speckleprop.IG_nucleoplasm=IG_nucleoplasm;
    speckleprop.IB_nucleoplasm=IB_nucleoplasm;
end

end