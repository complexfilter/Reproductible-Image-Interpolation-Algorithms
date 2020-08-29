

%%%%%%%  This function returns the block that we use for orthogonal block matching
%%%%%%%  pursuit. BNumber=28 types of blocks, including 20 orientations and length=6, 8,
%%%%%%%  12 are defined. The blocks are half-overlapped in the horizontal direction 
%%%%%%%  or vertical direction that is closer to the block orientation and
%%%%%%%  non-overlapped in the other direction. Blocks are shifted for
%%%%%%%  num_shift=12 times and the corresponding 12 block positions are
%%%%%%%  stored in BPosition. 
%%%%%%%
%%%%%%%  The blocks are stored in a semi-analytical way: for each type of
%%%%%%%  block, the geometry of only one "typical" block located near the
%%%%%%%  origin is stored. Other blocks of the same type can be obtained by
%%%%%%%  translating the typical block on a grid defined in BPosition. 
%%%%%%%
%%%%%%%  Input:
%%%%%%%  W----------- Width of the blocks. 
%%%%%%%  BNumber----- Number of block types (it should be 28).
%%%%%%%  N1,N2------- The image size, we need this while transforming the
%%%%%%%               row-column subscript to linear index for the points 
%%%%%%%               of a block.
%%%%%%%
%%%%%%%  Output:
%%%%%%%  BGeo-------- The block geometry. We note all points of a typical 
%%%%%%%               block with their linear index in the matrix, and this
%%%%%%%               typical block is located in a rectangle
%%%%%%%               [(1,1), (Support1,1)] * [(1,Support2),(Support1,Support2)]
%%%%%%%               BGeo is a bi-dimensional matrix size [L*W, BNumber].
%%%%%%%  BSize------- An array of size BNumber*1. BSize(k) is the number of points that  
%%%%%%%               the k-th tupe block has.
%%%%%%%  BPosition--- A matrix of size [N1*N2, BNumber, num_shift], with binary value.
%%%%%%%               If we can put the k-th type of block on
%%%%%%%               the n-th point in a image[N1,N2], and all this 
%%%%%%%               block's points are included strictly in this image,
%%%%%%%               then BPostion[n,k]=1. On the other hand, if any
%%%%%%%               point's subscript exceed the image's border,
%%%%%%%               BPosition[n,k] = 0.
%%%%%%%  Support1, Support2
%%%%%%%  ------------ the support of all the blocks, that means all the
%%%%%%%               blocks are contained in a rectangle Support1*Support2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BGeo,BSize,BPosition,Support1,Support2] = Define_Block_Dictionary(W,BNumber,N1,N2,num_shift)
BSize = zeros(BNumber,1);
% BPosition(n,k) stores the if the block of type k is defined on the
% coordinate n (BPosition(n,k) = 1) or not (BPosition(n,k) = 0).
% Compute BPosition of all the shifted grid. 
BPosition = false(N1*N2, BNumber, num_shift);
% S1(k) x S2(k) defines the minimum rectangular support that contains the
% block of type k positioned at the origin. 
S1 = zeros(BNumber,1);
S2 = zeros(BNumber,1);
L = 6;
Lmax = 12;
% This array stores the geometry of blocks positioned at the origin. 
BGeo = zeros(Lmax*W, BNumber);
% block of form ******** 
%               ******** 
k = 1;
for i = 1:L*W
   r = mod(i,W);
   if (r==0)
       r = W;
   end
   BGeo(i,k) =  N1*(i-r)/W+r;
end
S1(k) = W;
S2(k) = L; 
BSize(k) = L*W;


% block of form     *** 
%                ******
%             ******
%          ******
%          ***
%          
k = 2;
m = 3;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*(i-r)/W+floor((L-1)/m)-floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;


% block of form     **
%                 ****
%               ****
%             ****    
%             **
k = 3;
m = 2;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*((i-r)/W)+floor((L-1)/m)-floor((i-1)/(m*W));
    
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;

    
% block of form         *
%                      **
%                     **
%                    **
%                   **
%                  **
%                 **
%                **
%                *
k = 4;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = L+(N1-1)*(i-r)/W+r;
end
S1(k) = L+W;
S2(k) = L;
BSize(k) = L*W;


% block of form **
%               **
%              **
%              **
%             **
%             **
%            **
%            **
k = 5;
m = 2;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((L-1)/m)-floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;


% block of form **
%               **
%               **
%              **
%              **
%              **
%             **
%             **
%             **
%            **
%            **
%            **
k = 6;
m = 3;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((L-1)/m)-floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;

% block of form **
%               **
%               **
%               **
%               **
%               **
%               **
%               **
k = 7;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = N1*(i-r)/L+r;
end
S1(k) = L;
S2(k) = W;
BSize(k) = L*W;

% block of form **
%               **
%               **
%                **
%                **
%                **
%                  **
%                  **
%                  **
%                   **
%                   **
%                   **
k = 8;
m = 3;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;

% block of form **
%               **
%                **
%                **
%                 **
%                 **
%                  **
%                  **
k = 9;
m = 2;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*(floor((r-1)/m)+(i-r)/L);
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;


% block of form *
%               **
%                **
%                 **
%                  **
%                   **
%                    **
%                     **
%                      *
k = 10;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = (N1+1)*(i-r)/W+r;
   
end
S1(k) = L+W;
S2(k) = L;
BSize(k) = L*W;
 
% block of form **
%               ****
%                 ****
%                   ****
%                     **
k = 11;
m = 2;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = N1*(i-r)/W+floor((i-1)/(m*W))+r;

end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;

% block of form ***
%               ******
%                  ******
%                     ******
%                        ***
k = 12;
m = 3;
for i = 1:L*W
   r = mod(i,W);
   if (r==0)
       r =W;
   end
   BGeo(i,k) = r+N1*(i-r)/W+floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;

% block of form     **** 
%               ********
%               ****
k = 13;
m = 4;
L = 8;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*(i-r)/W+floor((L-1)/m)-floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;
% block of form **
%               **
%               **
%               **
%              **
%              **
%              **
%              **
k = 14;
m = 4;
L = 8;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r=L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((L-1)/m)-floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;
% block of form **
%               **
%               **
%               **
%                **
%                **
%                **
%                **
k = 15;
m = 4;
L = 8;
for i = 1:L*W
   r = mod(i,L);
   if (r==0)
       r = L;
   end
   BGeo(i,k) = r+N1*((i-r)/L+floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;
% block of form **** 
%               ********
%                   ****
k = 16;
m = 4;
L = 8;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*(i-r)/W+floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;
% block of form      ***** 
%               **********
%               *****
k = 17;
m = 6;
L = 12;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*(i-r)/W+floor((L-1)/m)-floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;
% block of form **
%               **
%               **
%               **
%              **
%              **
%              **
%              **
k = 18;
m = 6;
L = 12;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r=L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((L-1)/m)-floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;
% block of form **
%               **
%               **
%               **
%                **
%                **
%                **
%                **
k = 19;
m = 6;
L = 12;
for i = 1:L*W
   r = mod(i,L);
   if (r==0)
       r = L;
   end
   BGeo(i,k) = r+N1*((i-r)/L+floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;
% block of form **** 
%               ********
%                   ****
k = 20;
m = 6;
L = 12;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*(i-r)/W+floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;
% block of form     *** 
%                ******
%             ******
%          ******
%          ***         
k = 21;
m = 3;
L = 12;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*(i-r)/W+floor((L-1)/m)-floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;
% block of form **
%               **
%               **
%              **
%              **
%              **
%             **
%             **
%             **
%            **
%            **
%            **
k = 22;
m = 3;
L = 12;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((L-1)/m)-floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;
% block of form **
%               **
%               **
%                **
%                **
%                **
%                  **
%                  **
%                  **
%                   **
%                   **
%                   **
k = 23;
m = 3;
L = 12;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;
% block of form ***
%               ******
%                  ******
%                     ******
%                        ***
k = 24;
m = 3;
L = 12;
for i = 1:L*W
   r = mod(i,W);
   if (r==0)
       r =W;
   end
   BGeo(i,k) = r+N1*(i-r)/W+floor((i-1)/(m*W));
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;
% block of form     **
%                 ****
%               ****
%             ****    
%           **** 
%         **** 
%         **  
k = 25;
m = 2;
L = 12;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = r+N1*((i-r)/W)+floor((L-1)/m)-floor((i-1)/(m*W));
    
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;
% block of form **
%               **
%              **
%              **
%             **
%             **
%            **
%            **
%           **
%           **
%          **
%          **
k = 26;
m = 2;
L = 12;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*((i-r)/L+floor((L-1)/m)-floor((r-1)/m));
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;
% block of form **
%               **
%                **
%                **
%                 **
%                 **
%                  **
%                  **
%                   **
%                   **
%                    **
%                    **
k = 27;
m = 2;
L = 12;
for i = 1:L*W
    r = mod(i,L);
    if (r==0)
        r = L;
    end
    BGeo(i,k) = r+N1*(floor((r-1)/m)+(i-r)/L);
end
S1(k) = L;
S2(k) = W+floor((L-1)/m);
BSize(k) = L*W;

% block of form **
%               ****
%                 ****
%                   ****
%                     ****
%                       ****
%                       ****
k = 28;
m = 2;
L = 12;
for i = 1:L*W
    r = mod(i,W);
    if (r==0)
        r = W;
    end
    BGeo(i,k) = N1*(i-r)/W+floor((i-1)/(m*W))+r;
end
S1(k) = W+floor((L-1)/m);
S2(k) = L;
BSize(k) = L*W;

redund_factor_L = 2;
redund_factor_W = 1;


for k = 1:BNumber
    if ( ( k >= 1 ) && ( k <= 4 ) )
        delta1 = 2/redund_factor_W;
        delta2 = 6/redund_factor_L;
    elseif ( ( k >= 5 ) && ( k <= 9 ) )
        delta1 = 6/redund_factor_L;
        delta2 = 2/redund_factor_W;
    elseif ( ( k >= 10 ) && ( k <= 12 ) )
        delta1 = 2/redund_factor_W;
        delta2 = 6/redund_factor_L;
    elseif ( ( k == 13 ) || ( k == 16 ) ) 
        delta1 = 2/redund_factor_W;
        delta2 = 8/redund_factor_L;
    elseif ( ( k == 14 ) || ( k == 15 ) )
        delta1 = 8/redund_factor_L;
        delta2 = 2/redund_factor_W;
    elseif ( ( k == 17 ) || ( k == 20 ) )
        delta1 = 2/redund_factor_W;
        delta2 = 12/redund_factor_L;
    elseif ( ( k == 18 ) || ( k == 19 ) )
        delta1 = 12/redund_factor_L;
        delta2 = 2/redund_factor_W;
    end
    
    for kk = 1 : num_shift
        ind_shift1 = kk;
        ind_shift_k = mod(ind_shift1, delta1*delta2);
        if ( ind_shift_k == 0 )
            ind_shift_k = delta1*delta2;
        end
%         [shift1, shift2] = ind2sub([delta1, delta2], ind_shift_k);
        shift1 = mod(ind_shift_k, delta1);
        if ( shift1 == 0 )
            shift1 = delta1;
        end
        shift2 = ceil( ind_shift_k / delta1 );

        [YI, XI] = meshgrid(shift2 : delta2 : N2, shift1 : delta1 : N1);
%         ind = sub2ind([N1, N2], XI, YI);
        ind = (YI - 1) * N1 + XI;
        BPosition(ind, k, kk) = 1;
        
        % remove the border where the blocks go beyond
        for i1 = N1-S1(k)+2:N1
            BPosition(i1+N1*(0:N2-1), k, kk) = 0;
        end
        for i2 = N2-S2(k)+2:N2
            BPosition((1:N1)+(i2-1)*N1, k, kk) =0;
        end
    end
end

% Support1 x Support2 defines the minimum rectangular support that contain
% all types of blocks positioned at the origin. 
Support1 = max(S1(:));
Support2 = max(S2(:));
end