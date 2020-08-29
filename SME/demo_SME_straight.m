close all;clear;clc;
% Read high-resolution the image
% f = double(imread('pepper.png'));
addpath('/Users/complexfilter/Documents/paper code');
I_ini=edge_generate_straight_offset_0131...
   (25,129,129,257,257,4,200,0,0,2);% ground truth image. 
% I_ini=edge_generate_straight_offset_0131...
%    (45,129,129,257,257,4,200,0,0,2);% ground truth image. 
f=double(I_ini(1:512,1:512));
imH = f;
% Downsample the high-resolution by a factor of 2x2. 
imL = imH(1:2:end, 1:2:end);
% Super-resolution image zooming with mixing sparse estimators. 
% By defaut, no extrapolation.
option.extrapol = 0;
% Defaut parameter of orthogonal block matching pursuit.  
option.lambda = 0.6;
option.T = 0;
% The Lagragian parameter and the orthogonal block matching pursuit
% stopping threshold. The defaut values can be used for all the
% pictures. (See the paper for more details about these parameters.)
option.null = 0;
itp_method = 'spline';
BNumber = 28;
num_angles = 20;
[N1,N2] = size(imL);
flag_N1odd = 0;
flag_N2odd = 0;
Nh1 = 2 * N1 - 1;
Nh2 = 2 * N2 - 1;
% 1 scale wavelet transform
ScaleNumber = 1;
% Wavelet = cwt2(imL,ScaleNumber);
Wavelet = Wavelet_Transform_TI_1scale(imL);
% get the image low-frequency components
Wavelet_lf = zeros(N1, N2, 3*ScaleNumber+1);
Wavelet_lf(:,:,3*(ScaleNumber-1)+1) = Wavelet(:,:,3*(ScaleNumber-1)+1);
%% How do you computer the wavelet coefficients??

obtain_higher_response_concise
Wavelet_lf1 = zeros(N1, N2, 3*ScaleNumber+1);
Wavelet_lf1(:,:,3*(ScaleNumber-1)+1)=ddd1(2:2:end,2:2:end);

% 
% Wavelet_full = Wavelet_Transform_TI_1scale(imH);
% Wavelet_full_lf = zeros(N1, N2, 3*ScaleNumber+1);
% Wavelet_full(:,:,3*(ScaleNumber-1)+1) = Wavelet_full(:,:,3*(ScaleNumber-1)+1);
% 

im_LR_lf  = Wavelet_InverseTransform_TI_1scale(Wavelet_lf);
im_LR_lf1 = Wavelet_InverseTransform_TI_1scale(Wavelet_lf1);

% the low-resolution version's interpolated image by zero-padding higher
% frequency coefficients. 

Wavelet_hf = Wavelet;
Wavelet_hf(:,:,3*(ScaleNumber-1)+1) = 0;
im_SR_hf_all = zeros(Nh1, Nh2, BNumber);
[X,Y]   = meshgrid(1:0.5:N2, 1:0.5:N1);
[CX,CY] = meshgrid(1:N2, 1:N1);
delta1_max = 6;
delta2_max = 2;
num_shift  = delta1_max*delta2_max;
M_a = zeros(N1, N2, num_angles);% num_angles=20.
% block width
W = 2;
% Define the block dictionary.
[BGeo,BSize,BPosition_all,Support1,Support2] = Define_Block_Dictionary(W, BNumber, N1, N2, num_shift);
% Output:
% BPosition_all: 65536x28x12
% Support1: 12
% Support2: 12
% Support1 x Support2 defines the minimum rectangular support that contain
% all types of blocks positioned at the origin. 
% BSize: [ 12 12 12 12 12 12 12 12 12 12 12 12 16 16 16 16 24 24 24 24 24 24 24 24 24 24 24 24]'
% BSize(k) is the number of points that the k-th type block has.
% BGeo: 24x28  


% Input:
% W=2
% Bnumber:28(it means there exists 28 types of blocks in this paper)
% N1:256 (size of LR image)
% N2:256 (size of LR image)
% numshift:12 
h = waitbar(0,'Orthogonal block matching pursuit in progress: Please wait...');
for ind_shift = 1 : num_shift     % num_shift=12
%     ind_shift
    for j = 1:ScaleNumber
        if (j<ScaleNumber)
            direction1 = 3*j-2;
            direction2 = 3*j-1;
            direction3 = 3*j;
        else
            direction1 = 3*j-1;
            direction2 = 3*j;
            direction3 = 3*j+1;
        end
        Signal1 = Wavelet(:,:,direction1);
        Signal2 = Wavelet(:,:,direction2);
        Signal3 = Wavelet(:,:,direction3);        
        % Take the block dictionary grid. 
        BPosition = BPosition_all(:,:,ind_shift);
        % Compute the mixing coefficients with orthogonal block matching pursuit.
        [BListSize,BList,Bepsilon] = ComputeMixingCoefs_OBMP(Signal1,Signal2,Signal3, BGeo,BSize,BPosition,Support1,Support2,BNumber,W,option);
        % Read the mixing coefficients of each block type
        [AngleMapLR, EpsilonMapLR] = ReadMixingCoefs(BList,Bepsilon,BListSize,BGeo,BSize, N1, N2);   
        % Aggregate the mixing coefficients of each angle. 
        % (There are 20 angles and 28 block types, out of which 8 types of 
        % blocks (of various lengths) have the same angle
        % as other blocks. The mixing coefficients corresponding to the 
        % same angle are aggregated.)
        for k = 1 : BNumber % aggregate 28 types of blocks into 20 direction range.
            if     ( k == 21 )
                k1 = 2;
            elseif ( k == 22 )
                k1 = 6;
            elseif ( k == 23 )
                k1 = 8;
            elseif ( k == 24 )
                k1 = 12;
            elseif ( k == 25 )
                k1 = 3;
            elseif ( k == 26 )
                k1 = 5;
            elseif ( k == 27 )
                k1 = 9;
            elseif ( k == 28 )
                k1 = 11;
            else
                k1 = k;
            end
            M_a(:,:,k1) = M_a(:,:,k1) + EpsilonMapLR .* ( AngleMapLR == k );
%             if k==10
%             disp([k1 k EpsilonMapLR(116,116) AngleMapLR(116,116)  (M_a(116,116,1)) ]);
%             end
%             if k==2
%             disp([k1 k EpsilonMapLR(116,116) AngleMapLR(116,116)  (M_a(116,116,2)) ]);
%             end
        end  
    end   
    %waitbar(ind_shift/num_shift, h);
end
close(h)
% Average the mixing coefficients with respect to the shifted grids. 
M_a = M_a / num_shift;
%%% Adaptive interpolator computation. (cf Equation 32 in the paper)
% Linear interpolation on the residue.
Wavelet_hf1 = zeros(N1, N2, 3*ScaleNumber+1);
Wavelet_hf1(:,:,direction1) = Wavelet_hf(:,:,direction1) .* (1 - sum(M_a, 3));
Wavelet_hf1(:,:,direction2) = Wavelet_hf(:,:,direction2) .* (1 - sum(M_a, 3));
Wavelet_hf1(:,:,direction3) = Wavelet_hf(:,:,direction3) .* (1 - sum(M_a, 3));
imH_hf = Wavelet_InverseTransform_TI_1scale(Wavelet_hf1);
imH_linear_hf = interp2(CX,CY,imH_hf,X,Y,'spline');
h = waitbar(0, 'Directional interpolation is in process: Please wait...');
for k = 1 : 20
    if       ( k == 1 )
        m = 0;
        n = 1;
    elseif ( ( k == 2 ) || ( k == 21 ) )
        m = 1;
        n = 3;
    elseif ( ( k == 3 ) || ( k == 25 ) )
        m = 1;
        n = 2;
    elseif ( k == 4 )
        m = 1;
        n = 1;
    elseif ( ( k == 5 ) || ( k == 26 ) )
        m = 2;
        n = 1;
    elseif ( ( k == 6 ) || ( k == 22 ) )
        m = 3;
        n = 1;
    elseif ( k == 7 )
        m = 1;
        n = 0;
    elseif ( ( k == 8 ) || ( k == 23 ) )
        m = 3;
        n = -1;
    elseif ( ( k == 9 ) || ( k == 27 ) )
        m = 2;
        n = -1;
    elseif ( k == 10 )
        m = 1;
        n = -1;
    elseif ( ( k == 11 ) || ( k == 28 ) )
        m = 1;
        n = -2;
    elseif ( ( k == 12 ) || ( k == 24 ) )
        m = 1;
        n = -3;
    elseif ( k == 13 )
        m = 1;
        n = 4;
    elseif ( k == 14 )
        m = 4;
        n = 1;
    elseif ( k == 15 )
        m = 4;
        n = -1;
    elseif ( k == 16 )
        m = 1;
        n = -4;
    elseif ( k == 17 )
        m = 1;
        n = 6;
    elseif ( k == 18 )
        m = 6;
        n = 1;
    elseif ( k == 19 )
        m = 6;
        n = -1;
    elseif ( k == 20 )
        m = 1;
        n = -6;
    end        
    Wavelet_hf1(:,:,direction1) = Wavelet_hf(:,:,direction1) .* M_a(:,:,k);
    Wavelet_hf1(:,:,direction2) = Wavelet_hf(:,:,direction2) .* M_a(:,:,k);
    Wavelet_hf1(:,:,direction3) = Wavelet_hf(:,:,direction3) .* M_a(:,:,k);
    im_hf1 = Wavelet_InverseTransform_TI_1scale(Wavelet_hf1);    
    % Directional interpolation (cf the Appendix of the paper)
    im_SR_hf_all(:,:,k) = DirInterp2_3step(im_hf1, m, n, itp_method);
    waitbar(k/20, h);
end
close(h);
% Linear interpolation on the wavelet low-frequency. 
imH_linear_lf = interp2(CX,CY,im_LR_lf,X,Y, itp_method);
imH_linear_lf1 = interp2(CX,CY,im_LR_lf1,X,Y, itp_method);

% Combine the three interpolated components. (cf Equation 32 in the paper)
imH_SR = sum(im_SR_hf_all, 3) + imH_linear_hf + imH_linear_lf;
% linear interpolation
figure
imagesc(imL)
truesize
colormap('gray')
title('LR Inputimage');
figure
imagesc(imH)
truesize
colormap('gray')
title('HR ground truth')
PSNR = get_PSNR(imH(1:end-1, 1:end-1), imH_SR);
tmp = sprintf('SME supre-resolution %.2f dB', PSNR);
disp(tmp)
imH_H=zeros(512);
imH_H(1:511,1:511)=imH_SR(1:511,1:511);

imH_SR1 = sum(im_SR_hf_all, 3) + imH_linear_hf + imH_linear_lf1;
imH_SR1(1:2:end,1:2:end)=f(1:2:end,1:2:end);
clear PSNR;
PSNR(imH_SR(100:400,100:400),f(100:400,100:400))
PSNR(imH_SR1(100:400,100:400),f(100:400,100:400))

% Gain: 
PSNR(imH_SR1(100:400,100:400),f(100:400,100:400))-PSNR(imH_SR(100:400,100:400),f(100:400,100:400))
%% It seems to me that the Gain is not significant. 

q1=find(M_a(:,:,1)~=0);dd1=zeros(256);
q2=find(M_a(:,:,2)~=0);dd2=zeros(256);

