function mask=threshmask_withTHD_YI(image,blurradius,threshint_max)
blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
normlog=mat2gray(log(blur), [0 threshint_max]);
thresh=graythresh(normlog);
mask=imbinarize(normlog,'adaptive', 'Sensitivity', 0.6);
mask=imfill(mask,'holes');
end