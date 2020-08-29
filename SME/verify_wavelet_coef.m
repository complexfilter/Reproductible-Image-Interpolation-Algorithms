close all;clear;
load qmf.mat;
%% There exists some complication going on in the process of 
%% computing wavelet coefficients. 
N1=256;N2=256;
I_ini=edge_generate_straight_offset_0131...
   (30,129,129,257,257,4,200,0,0,2);% ground truth image. 
x=double(I_ini(1:512,1:512));
x=x(1:2:end,1:2:end);
f = double(imread('pepper.png'));
x = f(1:2:end,1:2:end);

    [n1, n2] = size(x);
    % 1-scale transform
    D = 1;

	tiwt = zeros((3*D+1)*n1, n2);
	lastx = (3*D*n1+1):(3*D*n1+n1); lasty = 1:n2;
	tiwt(lastx,lasty) = x;

    % l is no use
    l = 0;
    d = 0;
    s = x;

    wc00 = FWT2_PO_1scale_lantao(s,qmf);
    wc01 = FWT2_PO_1scale_lantao(CircularShift_lantao(s,0,1),qmf);
    wc10 = FWT2_PO_1scale_lantao(CircularShift_lantao(s,1,0),qmf);
    wc11 = FWT2_PO_1scale_lantao(CircularShift_lantao(s,1,1),qmf);
     
gg=cconv2(qmf(end:-1:1),qmf(end:-1:1),x,-4,-4);

figure;imagesc(circshift(gg(2:2:end,2:2:end),[-1 -1])-wc00(1:(n1/2),1:(n2/2)));
% note: circshift(gg(2:2:end,2:2:end),[-1 -1])=wc00(1:(n1/2),1:(n2/2)).
figure;imagesc(circshift(gg(1:2:end,1:2:end),[-1 -1])-wc11(1:(n1/2),1:(n2/2)));
% note: circshift(gg(1:2:end,1:2:end),[-1 -1])=wc11(1:(n1/2),1:(n2/2)).
figure;imagesc(circshift(gg(1:2:end,2:2:end),[-1 -1])-wc01(1:(n1/2),1:(n2/2)));
% note: circshift(gg(1:2:end,2:2:end),[-1 -1])=wc01(1:(n1/2),1:(n2/2)).
figure;imagesc(circshift(gg(2:2:end,1:2:end),[-1 -1])-wc10(1:(n1/2),1:(n2/2)));
% note: circshift(gg(2:2:end,1:2:end),[-1 -1])=wc10(1:(n1/2),1:(n2/2)).


% arrange the wavelet coefs to the usual setting
% CWTcoef = zeros(N1, N2);
% 
% CWTcoef(1 : 2 : end, 1 : 2 : end) = circshift(gg(1:2:end,1:2:end),[-1 -1]);
% CWTcoef(2 : 2 : end, 1 : 2 : end) = circshift(gg(2:2:end,1:2:end),[-1 -1]);
% CWTcoef(1 : 2 : end, 2 : 2 : end) = circshift(gg(1:2:end,2:2:end),[-1 -1]);
% CWTcoef(2 : 2 : end, 2 : 2 : end) = circshift(gg(2:2:end,2:2:end),[-1 -1]);

CWTcoef=cconv2(qmf(end:-1:1),qmf(end:-1:1),x,-6,-6);

load Wavelet_low_f.mat
figure;imagesc(CWTcoef-Wavelet_low_f);
norm(CWTcoef-Wavelet_low_f,'fro')
norm(circshift(gg,[-2 -2])-Wavelet_low_f,'fro')

