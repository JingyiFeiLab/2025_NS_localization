function output_img=autocontrast(img,r1,r2)
low_limit=r1;
up_limit=r2;
[m1 n1]=size(img);
img=double(img);
%--------------------calculation of vmin and vmax----------------------

    arr=sort(reshape(img,m1*n1,1));
    v_min=arr(ceil(low_limit*m1*n1));
    v_max=arr(ceil(up_limit*m1*n1));

%----------------------------------------------------------------------
% if r1==3
%     v_min=rgb2ntsc(v_min);
%     v_max=rgb2ntsc(v_max);
% end
%----------------------------------------------------------------------
output_img=(img-v_min)/(v_max-v_min);
output_img=uint16(65536*output_img);
end