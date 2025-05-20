function img_DoG = DoG(img,sigma1,sigma2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% sigma1<sigma2
im1=imgaussfilt(img,sigma1);
im2=imgaussfilt(img,sigma2);
img_DoG=im1-im2;
end