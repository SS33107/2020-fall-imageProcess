%%  Figure of s = T(r) 
r = [0:255];
sori = atan(double((r-128)/32));
s = ((sori-sori(1)) * 255) / (sori(256)-sori(1));
subplot(2,2,1);
plot(r,s,'b');
title('s = T(r) ');
%%  Table of transformation function 
table = [r;s];
T = array2table(table.','VariableNames',{'r','s'})
fig = uifigure;
uit = uitable(fig,'Data',T);
%writetable(T,'Table of transformation function .xlsx')
%%   output image 
I = imread('Bird feeding 3 low contrast.tif');
out = (atan((double(I)-128)/32)-sori(1))*(255/(sori(256)-sori(1)));
%1.3258
subplot(2,2,2);
imshow(out,[0 255]);
title('output image');
%%    original and output histograms
subplot(2,2,3);
histogram(I);
title('original histogram');
subplot(2,2,4);
histogram(out);
title('output histogram');