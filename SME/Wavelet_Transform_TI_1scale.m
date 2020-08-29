function CWTcoef = Wavelet_Transform_TI_1scale(f)
%
% 1-scale 2D TI wavelet transform (Daubechies wavelets with 4 vanishing
% moments). 
%
% Input: 
% f: an image of size NxN. (N must be 2^K, where K is an integer.)
%
% Output:
% CWTcoef: the wavelet coefficients NxNx4. CWTcoef(:,:,1:4) are
% respectively the low-frequency coefficients, the wavelet coefficients of
% horizontal, vertical and diagonal bands.
%
% Wrapper to Gabriel Peyre's perform_wavelet_transform, which itself
% wrappers to some Wavelab routines. 
%
% See also: Wavelet_InverseTransform_TI_1scale
% 
% Guoshen Yu, 2009

[N1, N2] = size(f);

NumScales = 1;

options.ti = 1;
% Jmin = log2(N1)-NumScales;

% Wavelet1 = perform_wavelet_transform(f, Jmin, +1, options);
Wavelet1 = perform_wavelet_transform_1scale(f, +1, options);
CWTcoef1 = zeros(N1, N2, 3*NumScales+1);

for k = 1 : 4    
    CWTcoef1(:,:,k) = Wavelet1((k-1)*N1+1 : k*N1, :);
end

% arrange the wavelet coefs to the usual setting
CWTcoef = zeros(N1, N2, 3*NumScales+1);

% % low-frequency
CWTcoef1_11 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 4);
CWTcoef1_21 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 4);
CWTcoef1_12 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 4);
CWTcoef1_22 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 4);
CWTcoef(1 : 2 : end, 1 : 2 : end, 1) = CWTcoef1_22;
CWTcoef(2 : 2 : end, 1 : 2 : end, 1) = CWTcoef1_12;
CWTcoef(1 : 2 : end, 2 : 2 : end, 1) = CWTcoef1_21;
CWTcoef(2 : 2 : end, 2 : 2 : end, 1) = CWTcoef1_11;

% horizontal band
CWTcoef1_11 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 2);
CWTcoef1_21 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 2);
CWTcoef1_12 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 2);
CWTcoef1_22 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 2);
CWTcoef(1 : 2 : end, 1 : 2 : end, 2) = CWTcoef1_22;
CWTcoef(2 : 2 : end, 1 : 2 : end, 2) = CWTcoef1_12;
CWTcoef(1 : 2 : end, 2 : 2 : end, 2) = CWTcoef1_21;
CWTcoef(2 : 2 : end, 2 : 2 : end, 2) = CWTcoef1_11;

% vertical band
CWTcoef1_11 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 1);
CWTcoef1_21 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 1);
CWTcoef1_12 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 1);
CWTcoef1_22 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 1);

CWTcoef(1 : 2 : end, 1 : 2 : end, 3) = CWTcoef1_22;
CWTcoef(2 : 2 : end, 1 : 2 : end, 3) = CWTcoef1_12;
CWTcoef(1 : 2 : end, 2 : 2 : end, 3) = CWTcoef1_21;
CWTcoef(2 : 2 : end, 2 : 2 : end, 3) = CWTcoef1_11;


% diagonal band
CWTcoef1_11 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 3);
CWTcoef1_21 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 3);
CWTcoef1_12 = CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 3);
CWTcoef1_22 = CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 3);

CWTcoef(1 : 2 : end, 1 : 2 : end, 4) = CWTcoef1_22;
CWTcoef(2 : 2 : end, 1 : 2 : end, 4) = CWTcoef1_12;
CWTcoef(1 : 2 : end, 2 : 2 : end, 4) = CWTcoef1_21;
CWTcoef(2 : 2 : end, 2 : 2 : end, 4) = CWTcoef1_11;
end
% % for k = 1 : 4
% %     figure
% %     imagesc(CWTcoef(:,:,k))
% %     colormap('gray')
% %     truesize
% % end
% % % % 
% % % % 
% % % % % f_rec = perform_wavelet_transform(Wavelet1, Jmin, -1, options);
% % % % 
