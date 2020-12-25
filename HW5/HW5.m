clc
clear 
close all
ori = imread('Car On Mountain Road.tif');
ori = im2double(ori);
%% Marr-Hildreth edge detection 
sigma = 4;
n = 25;
kernel_size = (25-1)/2;
[x, y] = meshgrid(-kernel_size : kernel_size, -kernel_size : kernel_size);
a = (x .^ 2 + y .^ 2 - 2 * sigma ^ 2) / sigma ^ 4;
b = exp( - (x .^ 2 + y .^ 2) / (2 * sigma ^ 2) );
LoG = a .* b;
LoG = LoG / sum(LoG(:));
figure
convResult = conv2(ori,LoG,'same');
imshow(convResult,[]);
[rr,cc]=size(convResult);
threshold4 = 0.04 * max(abs(convResult(:)));
zc0=zeros([rr,cc]);
zc4=zeros([rr,cc]);
for i=2:rr-1
    for j=2:cc-1
       if ((convResult(i,j+1)*convResult(i,j)<0) || (convResult(i,j)*convResult(i,j-1)<0)) 
          zc0(i,j)=1;
          if(abs(convResult(i,j+1)-convResult(i,j))>threshold4 || abs(convResult(i,j)-convResult(i,j-1))>threshold4)
              zc4(i,j)=1;
          end
       elseif((convResult(i+1,j)*convResult(i,j)<0) || (convResult(i,j)*convResult(i-1,j)<0)) 
          zc0(i,j)=1;
          if(abs(convResult(i+1,j)-convResult(i,j))>threshold4 || abs(convResult(i,j)-convResult(i-1,j))>threshold4)
              zc4(i,j)=1;
          end
       elseif((convResult(i,j)==0)&&(convResult(i,j+1)~=convResult(i,j-1))) 
           zc0(i,j)=1;
           if(abs(convResult(i,j+1)-convResult(i,j-1))>2*threshold4 )
              zc4(i,j)=1;
           end
         elseif((convResult(i,j)==0)&&(convResult(i+1,j)~=convResult(i-1,j))) 
           zc0(i,j)=1;
           if(abs(convResult(i+1,j)-convResult(i-1,j))>2*threshold4)
              zc4(i,j)=1;
          end
       end
    end
end
figure
imshow(zc0,[]);
figure
imshow(zc4,[]);
%% Edge linking by Hough transform
[H,T,R] = hough(zc4);
figure
imshow(H,[0 100],'XData',T,'YData',R,...
            'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
lines = houghlines(zc4,T,R,P,'FillGap',5,'MinLength',7);
[m,n]=size(zc4);
line = zeros(m,n);
line(:,:)=1;
figure, imshow(line), hold on
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end
figure, imshow(ori), hold on
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end