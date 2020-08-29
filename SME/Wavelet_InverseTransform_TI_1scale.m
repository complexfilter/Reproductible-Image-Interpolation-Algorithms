function f = Wavelet_InverseTransform_TI_1scale(CWTcoef)
%
% 1-scale 2D TI wavelet inverse transform (Daubechies wavelets with 4 vanishing
% moments). 
%
% Input: 
% CWTcoef: the wavelet coefficients NxNx4. CWTcoef(:,:,1:4) are
% respectively the low-frequency coefficients, the wavelet coefficients of
% horizontal, vertical and diagonal bands. It should be produced by
% Wavelet_Transform_TI_1scale
%
% Output:
% f: an image of size NxN. (N must be 2^K, where K is an integer.)
%
% Wrapper to Gabriel Peyre's perform_wavelet_transform, which itself
% wrappers to some Wavelab routines. 
%
% See also: Wavelet_Transform_TI_1scale
% 
% Guoshen Yu, 2009

N1 = size(CWTcoef, 1);
N2 = size(CWTcoef, 2);

NumScales = 1;

options.ti = 1;
% Jmin = log2(N1)-NumScales;

% arrange the wavelet coefs to the usual setting
CWTcoef1 = zeros(N1, N2, 3*NumScales+1);

% % low-frequency
CWTcoef1_22 = CWTcoef(1 : 2 : end, 1 : 2 : end, 1);
CWTcoef1_12 = CWTcoef(2 : 2 : end, 1 : 2 : end, 1);
CWTcoef1_21 = CWTcoef(1 : 2 : end, 2 : 2 : end, 1);
CWTcoef1_11 = CWTcoef(2 : 2 : end, 2 : 2 : end, 1);
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 4) = CWTcoef1_11;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 4) = CWTcoef1_21;
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 4) = CWTcoef1_12;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 4) = CWTcoef1_22;

% horizontal band
CWTcoef1_22 = CWTcoef(1 : 2 : end, 1 : 2 : end, 2);
CWTcoef1_12 = CWTcoef(2 : 2 : end, 1 : 2 : end, 2);
CWTcoef1_21 = CWTcoef(1 : 2 : end, 2 : 2 : end, 2);
CWTcoef1_11 = CWTcoef(2 : 2 : end, 2 : 2 : end, 2);
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 2) = CWTcoef1_11;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 2) = CWTcoef1_21;
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 2) = CWTcoef1_12;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 2) = CWTcoef1_22;

% vertical band
CWTcoef1_22 = CWTcoef(1 : 2 : end, 1 : 2 : end, 3);
CWTcoef1_12 = CWTcoef(2 : 2 : end, 1 : 2 : end, 3);
CWTcoef1_21 = CWTcoef(1 : 2 : end, 2 : 2 : end, 3);
CWTcoef1_11 = CWTcoef(2 : 2 : end, 2 : 2 : end, 3);
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 1) = CWTcoef1_11;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 1) = CWTcoef1_21;
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 1) = CWTcoef1_12;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 1) = CWTcoef1_22;


% diagonal band
CWTcoef1_22 = CWTcoef(1 : 2 : end, 1 : 2 : end, 4);
CWTcoef1_12 = CWTcoef(2 : 2 : end, 1 : 2 : end, 4);
CWTcoef1_21 = CWTcoef(1 : 2 : end, 2 : 2 : end, 4);
CWTcoef1_11 = CWTcoef(2 : 2 : end, 2 : 2 : end, 4);
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (1-1)*N2/2+1 : 1*N2/2, 3) = CWTcoef1_11;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (1-1)*N2/2+1 : 1*N2/2, 3) = CWTcoef1_21;
CWTcoef1((1-1)*N1/2+1 : 1*N1/2, (2-1)*N2/2+1 : 2*N2/2, 3) = CWTcoef1_12;
CWTcoef1((2-1)*N1/2+1 : 2*N1/2, (2-1)*N2/2+1 : 2*N2/2, 3) = CWTcoef1_22;

Wavelet1 = zeros(4*N1, N2);
for k = 1 : 4    
    Wavelet1((k-1)*N1+1 : k*N1, :) = CWTcoef1(:,:,k);
end

% f = perform_wavelet_transform(Wavelet1, Jmin, -1, options);
f = perform_wavelet_transform_1scale(Wavelet1, -1, options);
end
