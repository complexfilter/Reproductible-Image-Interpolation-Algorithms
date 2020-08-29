%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%%%%%%%%    Read the orientations and the mixing coefficients computed by
%%%%%%%%    ComputeMixingCoefs_OBMP and store them the corresponding places in
%%%%%%%%    the image. 
%%%%%%%%
%%%%%%%%    Input:
%%%%%%%%    BList,Bepsilon,BListSize: block information (orientations and mixing coefficients)
%%%%%%%%    calculated by ComputeMixingCoefs_OBMP.
%%%%%%%%    BGeo,BSize,: block dictionary parameters (see Define_Block_Dictionary).
%%%%%%%%    N1, N2: low-resolution image size. 
%%%%%%%%
%%%%%%%%    Output:
%%%%%%%%    AngleMapLR: an image that stores the block orientation. 
%%%%%%%%    EpsilonMapLR: an image that stores the block mixing
%%%%%%%%    coefficients. 
%%%%%%%%
%%%%%%%%    Reference:
%%%%%%%%    S.Mallat and G.Yu, Super-Resolution with Sparse Mixing Estimators,
%%%%%%%%    submitted to IEEE Trans. on Image Processing, 2009.
%%%%%%%%    Website: http://www.cmap.polytechnique.fr/~mallat/SME
%%%%%%%% 
%%%%%%%%    Contact: Guoshen Yu  yu@cmap.polytechnique.fr
%%%%%%%%
%%%%%%%%    Copyright 2009 Guoshen Yu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AngleMapLR, EpsilonMapLR] = ReadMixingCoefs(BList,Bepsilon,BListSize,BGeo,BSize, N1, N2)

AngleMapLR = zeros(N1,N2);
EpsilonMapLR = zeros(N1,N2);

for i = 0 : BListSize-1
    k = BList(BListSize-i,2);
    
    % read the angle and the mixing coefficients and place them in the LR image
    V_ind_LR = BList(BListSize-i,1)-1+BGeo(1:BSize(k),k);
    AngleMapLR(V_ind_LR) = k;
    epsilon = Bepsilon(BListSize-i,2);
    EpsilonMapLR(V_ind_LR) = epsilon;    
end
end
