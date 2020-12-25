%%  Plot of DFT magnitude in Log scale
a = imread('Bird 2.tif');
b = im2double(a);
[m,n] = size(b);
d = zeros(m,n);
for i = 1:m
    for j = 1:n
        d(i,j) = b(i,j).*(-1).^(i + j);
    end
end
e = fft2(d);
F2 = log(abs(e));
subplot(2,1,1);
imshow(F2,[]); colorbar;title('DFT magnitude in Log scale');

%%  Image constructed by DFT
z = zeros(m,n);
for i = 1:m
    for j = 1:n
        z(i,j) = sqrt((i-m/2).^2 + (j-n/2).^2);
    end
end
H = zeros(m,n);
K = zeros(m,n);
for i = 1:m
    for j = 1:n
        if z(i,j) <= 30  
            H(i,j) = 1;
            K(i,j) = 0;
        else
            H(i,j) = 0;
            K(i,j) = 1;
        end
    end
end
h1 = e.*H;
h2 = ifft2(h1);
h3 = zeros(m,n);
k1 = e.*K;
k2 = ifft2(k1);
k3 = zeros(m,n);
for i = 1:m
    for j = 1:n
        h3(i,j) = h2(i,j).*((-1).^(i+j));
        k3(i,j) = k2(i,j).*((-1).^(i+j));
    end
end
subplot(2,1,2);
imshow([b real(h3) real(k3)]);
title('input image                 radius < 30 pixels                      radius > 30 pixels ');
imwrite(real(h3),'h3.jpg');
imwrite(real(k3),'k3.jpg');
%%  Table of top 25 DFT frequencies
mag = zeros(m,n/2);
mag = (abs(e(:,1:n/2)));
magSort = sort(mag(:),'descend');
magList = zeros(25,3);
for k = 1:25
    for i = 1:m
       for j = 1:n/2
         if(mag(i,j)==magSort(k))
             magList(k,1) = e(i,j);
             magList(k,2) = i;
             magList(k,3) = j;
          end
        end
    end
end
T = array2table(magList,'VariableNames',{'abs(magnitude)','u','v'});
fig = uifigure;
uit = uitable(fig,'Data',T);