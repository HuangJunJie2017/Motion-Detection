function img_out=img_desample(img_src,sample_rate)
[w,h,~]=size(img_src);
img_out= img_src(sample_rate:sample_rate:w-1,sample_rate:sample_rate:h-1,:);