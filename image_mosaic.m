clc
clear all
close all
fprintf('Computer Vision, Assignment, 2010-2011\n');
fprintf('Title: Image Mosaic\n');
fprintf('Author: Dong Hao\n');
fprintf('Supervisor: Bogdan Matuszewski\n');

%% Calculate the matrix
[fname pname index]=uigetfile({'D:\DH research\Computer Vision\EL3105 CV 2010-2011\Assignment 2011\.jpg'},'Choose a moving image');
if index
    str=[pname fname];
    im1=imread(str);
end
[fname pname index]=uigetfile({'D:\DH research\Computer Vision\EL3105 CV 2010-2011\Assignment 2011\.jpg'},'Choose a static image');
if index
    str=[pname fname];
    im2=imread(str);
end

[r1 c1 d1]=size(im1);
[r2 c2 d2]=size(im2);
NumPoints=4;

image(im1),axis image; hold on;
[y1 x1]=ginput(NumPoints);
figure,image(im2),axis image; hold on;
[y2 x2]=ginput(NumPoints);

a=[];
for i=1:NumPoints % x2 y2 static
    a=[ a;
        x1(i) y1(i) 1 0 0 0 -x1(i)*x2(i) -y1(i)*x2(i) -x2(i);...
        0 0 0 x1(i) y1(i) 1 -x1(i)*y2(i) -y1(i)*y2(i) -y2(i)];
end
[U S V]=svd(a);
v=V(:,9);
m=reshape(v,3,3)';


%% mosaic  x=r y=c
c  = m*[ 1 1  r1 r1 ;
         1 c1 1  c1 ;
         1 1  1  1 ];
rmin=min((c(1,:)./c(3,:)));
cmin=min((c(2,:)./c(3,:)));
rmax=max((c(1,:)./c(3,:)));
cmax=max((c(2,:)./c(3,:)));
Sr=round(abs(min([rmin 1])));
Sc=round(abs(min([cmin 1])));
Srmax=round(max([rmax r2]));
Scmax=round(max([cmax c2]));

%%
% %mosaic without interpolation using source to target
% Im_out(Sr+1:Sr+r2,Sc+1:Sc+c2,:)=im2(1:r2,1:c2,:);
% for r=1:r1-1
%     for c=1:c1-1
%         k=m*[r c 1]';
%         x=round(k(1)/k(3));
%         y=round(k(2)/k(3));
%         Im_out(x+Sr+1,y+Sc+1,:)=im1(r,c,:);
%     end
% end
% figure,image(Im_out),axis image; hold on;
% title('mosaic without interpolation using source to target');

%% 
% % %mosaic without interpolation using target to source
% % clear Im_out
% % im1=im2uint8(im1);
% % im2=im2uint8(im2);
% % Im_out(Sr+1:Sr+r2,Sc+1:Sc+c2,:)=im2(1:r2,1:c2,:);
% % for r=1:Srmax+Sr
% %     for c=1:Scmax+Sc
% %         k=inv(m)*[r-Sr c-Sc 1]';
% %         x=round(k(1)/k(3));
% %         y=round(k(2)/k(3));
% %         if (x>1 && x<r1) && (y>1 && y<c1) % error here
% %             Im_out(r,c,:)=im1(x,y,:);
% %         end
% %     end
% % end
% % figure,image(Im_out),axis image; hold on;
% % title('mosaic without interpolation using target to source');


%% mosaic with interpolation
clear Im_out
im1=im2double(im1);
im2=im2double(im2);
[cm,rm]=meshgrid(1-Sc:Scmax,1-Sr:Srmax);
[r,c]=size(rm);        %equals size(Yp)

K=inv(m)*[rm(:)'; cm(:)'; ones(1,r*c)] ; %warp
ro=reshape( K(1,:)./K(3,:),r,c);
co=reshape( K(2,:)./K(3,:),r,c);

[rm,cm]=meshgrid(1:c1,1:r1);

Im_out(:,:,1)=interp2(rm,cm,im1(:,:,1),co,ro, 'linear'); % red
Im_out(:,:,2)=interp2(rm,cm,im1(:,:,2),co,ro, 'linear'); % green
Im_out(:,:,3)=interp2(rm,cm,im1(:,:,3),co,ro, 'linear'); % blue

Im_out(Sr+1:Sr+r2,Sc+1:Sc+c2,:)=im2(1:r2,1:c2,:);
figure,image(Im_out),axis image; hold on;
title('mosaic with interpolation using target to source');