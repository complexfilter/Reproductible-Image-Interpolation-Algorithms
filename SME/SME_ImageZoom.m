function [imH_SR, imH_linear] = SME_ImageZoom(imL, option)
%
% Image zooming by a factor of 2x2 with sparse mixing estimators. 
%
% Input: 
% imL: gray-level low-resolution image. It should be a 2D matrix. 
%
% option (optional): it specifies two paramters in the mixing coefficient computation 
% with orthogonal block matching pursuit, as well as an extrapolation flag.
%   option.lambda defines the Lagragian parameter (defaut value 0.6). 
%
%   option.T is a threshold that defines the matching pursuit stopping criterion. 
%   It should be typically determined by the noise energy (defaut value 0). For image
%   zooming (when there is no noise), option.T=0 is appropriate.
% 
%   option.extrapol (defaut value 0): if option.extrapol=1, extrapolate the last row and
%   last column by line/column duplication. Let N1xN2 be the size of the
%   low-resolution image. If option.extrapol=1, the zoomed image with be of
%   size (2xN1)x(2xN2). Otherwise the zoomed image with be of size (2xN1-1)x(2xN2-1).
%   
% Output: 
% imH_SR: super-resolution image estimation obtained with sparse mixing estimators. 
% 
% imH_linear: cubic spline linear interpolation (for comparison).
%
% See also:
% Define_Block_Dictionary, ComputeMixingCoefs_OBMP, DirInterp2_3step
%
% Reference
%          
% S.Mallat and G.Yu, Super-Resolution with Sparse Mixing Estimators,
% submitted to IEEE Trans. on Image Processing, 2009.
% Website: http://www.cmap.polytechnique.fr/~mallat/SME
% 
% Contact: Guoshen Yu  yu@cmap.polytechnique.fr
%
% Copyright 2009 Guoshen Yu
% 
if ( ndims(imL) ~= 2 )
    error('The input image should be in gray-level and be a 2D matrix.');
elseif ( ( size(imL, 1) == 1 ) || ( size(imL, 2) == 1 ) )
    error('The input image should be in gray-level and be a 2D matrix.');
end

% The Lagragian parameter and the orthogonal block matching pursuit
% stopping threshold. The defaut values can be used for all the
% pictures. (See the paper for more details about these parameters.)
option.null = 0;
if isfield(option, 'lambda')
    if ( ~isnumeric(option.lambda) )
        error('option.lambda should be a number that specifies the Lagrangian value. By defaut option.lambda = 0.6.');
    end
else
    % by defaut, lambda = 0.6
    option.lambda = 0.6;
end

if isfield(option, 'T')
    if ( ~isnumeric(option.T) )
        error('option.T should be a number that specifies the thresholding value. By defaut option.T = 0.');
    end
else
    % by defaut, threshold = 0
    option.T = 0;
end

if isfield(option, 'extrapol')
    if ( ~isnumeric(option.extrapol) )
        error('option.extrapol should be 0 (no extrapolation) or 1 (extrapolating the last row and last column). By defaut option.extrapol = 0.');
    end
else
    % by defaut, threshold = 0
    option.extrapol = 0;
end

itp_method = 'spline';
%
BNumber = 28;
num_angles = 20;

[N1,N2] = size(imL);

% If the image size is odd, add a row/column to make it even. (The 1-scale wavelet transform
% requres that the image size is even.)
flag_N1odd = 0;
flag_N2odd = 0;
if ( mod(N1, 2) == 1 )
    tmp = zeros(N1+1, N2);
    tmp(1:N1, :) = imL;
    tmp(end, :) = imL(end, :);
    imL = tmp;
    N1 = N1 + 1;
    flag_N1odd = 1;
end
if ( mod(N2, 2) == 1 )
    tmp = zeros(N1, N2+1);
    tmp(:, 1:N2) = imL;
    tmp(:, end) = imL(:, end);
    imL = tmp;
    N2 = N2 + 1;
    flag_N2odd = 1;
end

Nh1 = 2 * N1 - 1;
Nh2 = 2 * N2 - 1;
% 1 scale wavelet transform
ScaleNumber = 1;
% Wavelet = cwt2(imL,ScaleNumber);
Wavelet = Wavelet_Transform_TI_1scale(imL);

% get the image low-frequency components
Wavelet_lf = zeros(N1, N2, 3*ScaleNumber+1);
Wavelet_lf(:,:,3*(ScaleNumber-1)+1) = Wavelet(:,:,3*(ScaleNumber-1)+1);
% im_LR_lf = cwt2i(Wavelet_lf);
im_LR_lf = Wavelet_InverseTransform_TI_1scale(Wavelet_lf);
% get the image high-frequency components
Wavelet_hf = Wavelet;
Wavelet_hf(:,:,3*(ScaleNumber-1)+1) = 0;
im_SR_hf_all = zeros(Nh1, Nh2, BNumber);
[X,Y] = meshgrid(1:0.5:N2, 1:0.5:N1);
[CX,CY] = meshgrid(1:N2, 1:N1);
delta1_max = 6;
delta2_max = 2;
num_shift = delta1_max*delta2_max;
M_a = zeros(N1, N2, num_angles);
% block width
W = 2;
% Define the block dictionary.
[BGeo,BSize,BPosition_all,Support1,Support2] = Define_Block_Dictionary(W, BNumber, N1, N2, num_shift);
% disp('Orthogonal block matching pursuit...')
% tic;
h = waitbar(0,'Orthogonal block matching pursuit in progress: Please wait...');
for ind_shift = 1 : num_shift    
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

        % Compute the mixing coefficients with orthogonal block matching
        % pursuit.
        [BListSize,BList,Bepsilon] = ComputeMixingCoefs_OBMP(Signal1,Signal2,Signal3, BGeo,BSize,BPosition,Support1,Support2,BNumber,W,option);

        % Read the mixing coefficients of each block type
        [AngleMapLR, EpsilonMapLR] = ReadMixingCoefs(BList,Bepsilon,BListSize,BGeo,BSize, N1, N2);
   
        % Aggregate the mixing coefficients of each angle. (There are 20 angles and 28 block types, out of which 8 types of blocks (of various lengths) have the same angle
        % as other blocks. The mixing coefficients corresponding to the same angle are aggregated.)
        for k = 1 : BNumber
            if ( k == 21 )
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
        end  
    end
    
    waitbar(ind_shift/num_shift, h);

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
% imH_hf = cwt2i(Wavelet_hf1);
imH_hf = Wavelet_InverseTransform_TI_1scale(Wavelet_hf1);
imH_linear_hf = interp2(CX,CY,imH_hf,X,Y,'spline');
% disp('Dir Interpolation')
% tic
% Directional interpolation on the directional components of different
% angles. 
h = waitbar(0, 'Directional interpolation is in process: Please wait...');
for k = 1 : 20
    if ( k == 1 )
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
%     im_hf1 = cwt2i(Wavelet_hf1);
     im_hf1 = Wavelet_InverseTransform_TI_1scale(Wavelet_hf1); 
    

    % Directional interpolation (cf the Appendix of the paper)
%     im_SR_hf_all(:,:,k) = DirInterp2_v3(im_hf1, m, n, itp_method);
    im_SR_hf_all(:,:,k) = DirInterp2_3step(im_hf1, m, n, itp_method);
    waitbar(k/20, h);
end
close(h)

% Linear interpolation on the wavelet low-frequency. 
imH_linear_lf = interp2(CX,CY,im_LR_lf,X,Y, itp_method);

% Combine the three interpolated components. (cf Equation 32 in the paper)
imH_SR = sum(im_SR_hf_all, 3) + imH_linear_hf + imH_linear_lf;

% linear interpolation
imH_linear = interp2(CX,CY,imL,X,Y, 'spline');

if ( flag_N1odd == 1 )
    imH_SR = imH_SR(1:end-2, :);
    imH_linear = imH_linear(1:end-2, :);
    Nh1 = Nh1 - 2;
end
if ( flag_N2odd == 1 )
    imH_SR = imH_SR(:, 1:end-2);
    imH_linear = imH_linear(:, 1:end-2);
    Nh2 = Nh2 - 2;
end

if ( option.extrapol == 1 )
    tmp = zeros(Nh1+1, Nh2+1);
    tmp(1:Nh1, 1:Nh2) = imH_SR;
    tmp(end, 1:Nh2) = imH_SR(end, 1:Nh2);
    tmp(1:Nh1, end) = imH_SR(1:Nh1, end);
    tmp(end, end) = imH_SR(end, end);
    imH_SR = tmp;
    
    tmp = zeros(Nh1+1, Nh2+1);
    tmp(1:Nh1, 1:Nh2) = imH_linear;
    tmp(end, 1:Nh2) = imH_linear(end, 1:Nh2);
    tmp(1:Nh1, end) = imH_linear(1:Nh1, end);
    tmp(end, end) = imH_linear(end, end);
    imH_linear = tmp;
end
end