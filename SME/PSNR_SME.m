close all;clear;clc;
%% The code below tests how the complex-valued curvelets help 
%% the interpolation.
addpath('/Users/complexfilter/Dropbox/Dropbox/code/Directional Processing');
addpath('/Users/complexfilter/Documents/paper code');
addpath('/Users/complexfilter/Dropbox/Dropbox/code/Modification/Label_Modification/');

grid_ini
grid_ini_B_0525
filter_ini_band_all_A
Forward_Filter_1x_2x_1210
sharpness=1;
PSNR_all=zeros(180,1);
indd=1;
for ang=-180:-1
xx=double(edge_generate_straight_offset_0131(ang,256,256,512,512,sharpness,100,50,0));
xx(xx>=255)=255;
xx(xx<0)=0;
I_observed=xx(1:2:end,1:2:end);
temp=abs(180+ang-90);
if abs(temp)>=45
y_b1=20.5:1:230.5;
if abs(temp)==90
    x_b1=128.5*ones(size(y_b1));
else
    x_b1=(y_b1-128.5)*tan(ang*pi/180)+128.5;
end
else
x_b1=20.5:1:230.5;
if abs(temp)==0
    y_b1=128.5*ones(size(x_b1));
else
    y_b1=(x_b1-128.5)/tan(ang*pi/180)+128.5;
end    
end
r=50;
orientation_field=zeros(256);
[~,region_ind_bb]=line2region(x_b1-0.5,y_b1-0.5,r,256,256);
region_ind_bb=round(region_ind_bb);
orientation_field(region_ind_bb)=180+ang;

r=50;
[~,region_ind_bb]=line2region(x_b1-0.5,y_b1-0.5,r,256,256);
RR=zeros(256);region_ind_bb=round(region_ind_bb);
RR(region_ind_bb)=1;
r=30;
[~,region_ind_bb]=line2region(x_b1-0.5,y_b1-0.5,r,256,256);
PP=zeros(256);region_ind_bb=round(region_ind_bb);
PP(region_ind_bb)=1;
r=20;
[~,region_ind_bb]=line2region(x_b1-0.5,y_b1-0.5,r,256,256);
QQ=zeros(256);region_ind_bb=round(region_ind_bb);
QQ(region_ind_bb)=1;
tic;
RECOS=sme(I_observed);
toc;
r=4;
[~,region_ind_bb]=line2region(2*x_b1-1,2*y_b1-1,r,512,512);
QQ_2x=zeros(512);
QQ_2x(region_ind_bb)=1;
mask=zeros(512);
mask(51:462,51:462)=1;
QQ_2x=QQ_2x.*mask;
PSNR_all(indd)=PSNR(RECOS(QQ_2x==1),xx(QQ_2x==1));
disp([indd-1 PSNR_all(indd)]);
indd=indd+1;
end

PSNR_all1=zeros(180,1);
PSNR_all1(1:46)=PSNR_all(1:46);
PSNR_all1(47:91)=PSNR_all(45:-1:1);
PSNR_all1(92:180)=PSNR_all1(90:-1:2);
