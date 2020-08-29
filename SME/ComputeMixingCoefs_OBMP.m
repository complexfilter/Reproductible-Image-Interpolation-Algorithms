%%%%%%%%    Compute the mixing coefficients with orthogonal block matching
%%%%%%%%    pursuit. 
%%%%%%%%
%%%%%%%%    Input:
%%%%%%%%    Signal1,Signal2,Signal3: the first-scale wavelet coefficients (of
%%%%%%%%    horizontal, vertical and diagonal bands) of the low-resolution
%%%%%%%%    image.
%%%%%%%%    BGeo,BSize,BPosition,Support1,Support2,W: block dictionary
%%%%%%%%    parameters (see Define_Block_Dictionary).
%%%%%%%%    option: (optional). option.lambda defines the Lagragian
%%%%%%%%    parameter (defaut value 0.6). option.T is a threshold that
%%%%%%%%    defines the matching pursuit stopping criterion. It should be typically
%%%%%%%%    determined by the noise energy (defaut value 0). For image
%%%%%%%%    zooming (when there is no noise), option.T=0 is appropriate.
%%%%%%%%
%%%%%%%%    Output:
%%%%%%%%    BListSize: number of blocks calculated. 
%%%%%%%%    BList: a matrix of size [BListSize, 2]. Each row stores the
%%%%%%%%    coordinate and the type of one block. 
%%%%%%%%    Bepsilon: a matrix of size [BListSize, 2]. Each row stores the
%%%%%%%%    coordinate and the mixing coefficient of one block (see Eq 26).
function [BListSize,BList,Bepsilon] = ComputeMixingCoefs_OBMP(Signal1,Signal2,Signal3, BGeo,BSize,BPosition,Support1,Support2,BNumber,W,option)

[N1,N2] = size(Signal1);
N = N1*N2;

% Initialisation
BListSize = 0;
BList = zeros(N,2);
Bepsilon = zeros(N,2);

Energy = repmat(0, [N,BNumber]);

Epsilon = repmat(0, [N,BNumber]);

% disp('Energy initiliazation...')


% alpha = -0.4;
% lambda = 1 + alpha;

option.null = 0;
if isfield(option, 'lambda')
    if ( ~isnumeric(option.lambda) )
        error('option.lambda should be a number that specifies the Lagrangian value. By defaut option.lambda = 0.6.');
    end
    lambda = option.lambda;
else
    % by defaut, lambda = 0.6
    lambda = 0.6;
end

if isfield(option, 'T')
    if ( ~isnumeric(option.T) )
        error('option.T should be a number that specifies the thresholding value. By defaut option.T = 0.');
    end
    Th = option.T;
else
    % by defaut, threshold = 0
    Th = 0;
end

% lambda = 0.6;
% Th = 0;

% loop over all types of blocks
for k = 1:BNumber
    % loop over all positions of a given type of block
    for n = 1:N
        % the block is defined on the position n if BPosition(n,k)=1
        if (BPosition(n,k))
            % take the wavelet coefficients in the block in the 3 wavelet
            % orientations
            BData1 = Signal1(n-1+BGeo(1:BSize(k),k));
            BData2 = Signal2(n-1+BGeo(1:BSize(k),k));
            BData3 = Signal3(n-1+BGeo(1:BSize(k),k));

            % Coefficient energy in the block    
            V_en_all = sum(BData1.^2) + sum(BData2.^2) + sum(BData3.^2);
            
            % Compute and sum the average energy on each block row
            V_en_avg = 0;
            for w = 1 : W
                if ( ( k == 1 ) || ( k == 3 ) || ( k == 4 ) || ( k == 10 ) || ( k == 11 ) )
                    V_en_avg = V_en_avg + ( (BData1(w)+BData1(W+w)+BData1(2*W+w)+ ...
                        BData1(3*W+w)+BData1(4*W+w)+BData1(5*W+w))^2/6 ) ...
                        + ( (BData2(w)+BData2(W+w)+BData2(2*W+w)+ ...
                        BData2(3*W+w)+BData2(4*W+w)+BData2(5*W+w))^2/6 ) ...
                        + ( (BData3(w)+BData3(W+w)+BData3(2*W+w)+ ...
                        BData3(3*W+w)+BData3(4*W+w)+BData3(5*W+w))^2/6 );

                elseif ( ( k == 2 ) || ( k == 12 ) )
                    V_en_avg = V_en_avg + ( (BData1(w)+BData1(W+w)+BData1(3*W+w)+ ...
                        BData1(2*W+w)+BData1(4*W+w)+BData1(5*W+w))^2/6 ) ...
                        + ( (BData2(w)+BData2(W+w)+BData2(3*W+w)+ ...
                        BData2(2*W+w)+BData2(4*W+w)+BData2(5*W+w))^2/6 ) ...
                        + ( (BData3(w)+BData3(W+w)+BData3(3*W+w)+ ...
                        BData3(2*W+w)+BData3(4*W+w)+BData3(5*W+w))^2/6 );

                elseif ( ( k == 5 ) || ( k == 7 ) || ( k == 9 ) )
                    L = 6;
                    V_en_avg = V_en_avg + ( (BData1((w-1)*L+1)+BData1((w-1)*L+2)+BData1((w-1)*L+3)+ ...
                        BData1((w-1)*L+4)+BData1((w-1)*L+5)+BData1((w-1)*L+6))^2/6 ) ...
                        + ( (BData2((w-1)*L+1)+BData2((w-1)*L+2)+BData2((w-1)*L+3)+ ...
                        BData2((w-1)*L+4)+BData2((w-1)*L+5)+BData2((w-1)*L+6))^2/6 ) ...
                        + ( (BData3((w-1)*L+1)+BData3((w-1)*L+2)+BData3((w-1)*L+3)+ ...
                        BData3((w-1)*L+4)+BData3((w-1)*L+5)+BData3((w-1)*L+6))^2/6 );

                elseif ( ( k == 6 ) || ( k == 8 ) )
                    L = 6;
                    V_en_avg = V_en_avg + ( (BData1((w-1)*L+1)+BData1((w-1)*L+2)+BData1((w-1)*L+4)+ ...
                        BData1((w-1)*L+3)+BData1((w-1)*L+5)+BData1((w-1)*L+6))^2/6 ) ...
                        + ( (BData2((w-1)*L+1)+BData2((w-1)*L+2)+BData2((w-1)*L+4)+ ...
                        BData2((w-1)*L+3)+BData2((w-1)*L+5)+BData2((w-1)*L+6))^2/6 ) ...
                        + ( (BData3((w-1)*L+1)+BData3((w-1)*L+2)+BData3((w-1)*L+4)+ ...
                        BData3((w-1)*L+3)+BData3((w-1)*L+5)+BData3((w-1)*L+6))^2/6 );

                elseif ( ( k == 13 ) || ( k == 16 ) )
                    V_en_avg = V_en_avg + ( (BData1(w)+BData1(W+w)+BData1(4*W+w)+BData1(5*W+w)+ ...
                        BData1(2*W+w)+BData1(3*W+w)+BData1(6*W+w)+BData1(7*W+w))^2/8 ) ...
                        + ( (BData2(w)+BData2(W+w)+BData2(4*W+w)+BData2(5*W+w)+ ...
                        BData2(2*W+w)+BData2(3*W+w)+BData2(6*W+w)+BData2(7*W+w))^2/8 ) ...
                        + ( (BData3(w)+BData3(W+w)+BData3(4*W+w)+BData3(5*W+w)+ ...
                        BData3(2*W+w)+BData3(3*W+w)+BData3(6*W+w)+BData3(7*W+w))^2/8 );

                elseif ( ( k == 14 ) || ( k == 15 ) )
                    L = 8;
                    V_en_avg = V_en_avg + ( (BData1((w-1)*L+1)+BData1((w-1)*L+2)+BData1((w-1)*L+5)+BData1((w-1)*L+6)+ ...
                        BData1((w-1)*L+3)+BData1((w-1)*L+4)+BData1((w-1)*L+7)+BData1((w-1)*L+8))^2/8 ) ...
                        + ( (BData2((w-1)*L+1)+BData2((w-1)*L+2)+BData2((w-1)*L+5)+BData2((w-1)*L+6)+ ...
                        BData2((w-1)*L+3)+BData2((w-1)*L+4)+BData2((w-1)*L+7)+BData2((w-1)*L+8))^2/8 ) ...
                        + ( (BData3((w-1)*L+1)+BData3((w-1)*L+2)+BData3((w-1)*L+5)+BData3((w-1)*L+6)+ ...
                        BData3((w-1)*L+3)+BData3((w-1)*L+4)+BData3((w-1)*L+7)+BData3((w-1)*L+8))^2/8 );

                elseif ( ( k == 17 ) || ( k == 20 ) )
                    V_en_avg = V_en_avg + ( (BData1(w)+BData1(W+w)+BData1(2*W+w)+BData1(6*W+w)+BData1(7*W+w)+BData1(8*W+w)+ ...
                        BData1(3*W+w)+BData1(4*W+w)+BData1(5*W+w)+BData1(9*W+w)+BData1(10*W+w)+BData1(11*W+w))^2/12) ...
                        + ( (BData2(w)+BData2(W+w)+BData2(2*W+w)+BData2(6*W+w)+BData2(7*W+w)+BData2(8*W+w)+ ...
                        BData2(3*W+w)+BData2(4*W+w)+BData2(5*W+w)+BData2(9*W+w)+BData2(10*W+w)+BData2(11*W+w))^2/12) ...
                        + ( (BData3(w)+BData3(W+w)+BData3(2*W+w)+BData3(6*W+w)+BData3(7*W+w)+BData3(8*W+w)+ ...
                        BData3(3*W+w)+BData3(4*W+w)+BData3(5*W+w)+BData3(9*W+w)+BData3(10*W+w)+BData3(11*W+w))^2/12);

                elseif ( ( k == 18 ) || ( k == 19 ) )
                    L = 12;
                    V_en_avg = V_en_avg + ( (BData1((w-1)*L+1)+BData1((w-1)*L+2)+BData1((w-1)*L+3)+BData1((w-1)*L+7)+BData1((w-1)*L+8)+BData1((w-1)*L+9)+  ...
                        BData1((w-1)*L+4)+BData1((w-1)*L+5)+BData1((w-1)*L+6)+BData1((w-1)*L+10)+BData1((w-1)*L+11)+BData1((w-1)*L+12) )^2/12 ) ...
                        + ( (BData2((w-1)*L+1)+BData2((w-1)*L+2)+BData2((w-1)*L+3)+BData2((w-1)*L+7)+BData2((w-1)*L+8)+BData2((w-1)*L+9)+  ...
                        BData2((w-1)*L+4)+BData2((w-1)*L+5)+BData2((w-1)*L+6)+BData2((w-1)*L+10)+BData2((w-1)*L+11)+BData2((w-1)*L+12) )^2/12 ) ...
                        + ( (BData3((w-1)*L+1)+BData3((w-1)*L+2)+BData3((w-1)*L+3)+BData3((w-1)*L+7)+BData3((w-1)*L+8)+BData3((w-1)*L+9)+  ...
                        BData3((w-1)*L+4)+BData3((w-1)*L+5)+BData3((w-1)*L+6)+BData3((w-1)*L+10)+BData3((w-1)*L+11)+BData3((w-1)*L+12) )^2/12 );

                elseif ( ( k == 21 ) || ( k == 24 ) || ( k == 25) || ( k == 28 ) )
                    V_en_avg = V_en_avg + ( (BData1(w)+BData1(W+w)+BData1(2*W+w)+BData1(6*W+w)+BData1(7*W+w)+BData1(8*W+w)+ ...
                        BData1(3*W+w)+BData1(4*W+w)+BData1(5*W+w)+BData1(9*W+w)+BData1(10*W+w)+BData1(11*W+w))^2/12) ...
                        + ( (BData2(w)+BData2(W+w)+BData2(2*W+w)+BData2(6*W+w)+BData2(7*W+w)+BData2(8*W+w)+ ...
                        BData2(3*W+w)+BData2(4*W+w)+BData2(5*W+w)+BData2(9*W+w)+BData2(10*W+w)+BData2(11*W+w))^2/12) ...
                        + ( (BData3(w)+BData3(W+w)+BData3(2*W+w)+BData3(6*W+w)+BData3(7*W+w)+BData3(8*W+w)+ ...
                        BData3(3*W+w)+BData3(4*W+w)+BData3(5*W+w)+BData3(9*W+w)+BData3(10*W+w)+BData3(11*W+w))^2/12);

                elseif ( ( k == 22 ) || ( k == 23 ) || ( k == 26 ) || ( k == 27 ) )
                    L = 12;
                    V_en_avg = V_en_avg + ( (BData1((w-1)*L+1)+BData1((w-1)*L+2)+BData1((w-1)*L+3)+BData1((w-1)*L+7)+BData1((w-1)*L+8)+BData1((w-1)*L+9)+  ...
                        BData1((w-1)*L+4)+BData1((w-1)*L+5)+BData1((w-1)*L+6)+BData1((w-1)*L+10)+BData1((w-1)*L+11)+BData1((w-1)*L+12) )^2/12 ) ...
                        + ( (BData2((w-1)*L+1)+BData2((w-1)*L+2)+BData2((w-1)*L+3)+BData2((w-1)*L+7)+BData2((w-1)*L+8)+BData2((w-1)*L+9)+  ...
                        BData2((w-1)*L+4)+BData2((w-1)*L+5)+BData2((w-1)*L+6)+BData2((w-1)*L+10)+BData2((w-1)*L+11)+BData2((w-1)*L+12) )^2/12 ) ...
                        + ( (BData3((w-1)*L+1)+BData3((w-1)*L+2)+BData3((w-1)*L+3)+BData3((w-1)*L+7)+BData3((w-1)*L+8)+BData3((w-1)*L+9)+  ...
                        BData3((w-1)*L+4)+BData3((w-1)*L+5)+BData3((w-1)*L+6)+BData3((w-1)*L+10)+BData3((w-1)*L+11)+BData3((w-1)*L+12) )^2/12 );

                end
            end
         
            Epsilon(n,k) = max( (1 - lambda*(V_en_all - V_en_avg)/V_en_all), 0 );
            Energy(n,k) = V_en_all * Epsilon(n,k)^2;

        end
    end

end
% ttt1=reshape(Energy(:,10),256,256);
% ttt=reshape(Epsilon(:,10),256,256);
% [ttt(116,116) ttt1(116,116)]
% if ttt1(116,116)>0.01
%     disp('stop');
% end
   
BlockChosen = zeros(N, 1);

% Sort the block energy
[Energy_sorted, IDX_sorted] = sort(Energy(:), 'descend');

NxBNumber = N * BNumber;


kk = 1;
while ( kk <= NxBNumber )
    % find the next energy (before the second last) that has not been set to zero
    while ( ( kk < NxBNumber ) && ( Energy(IDX_sorted(kk)) == 0 ) )
        kk = kk + 1;
    end
%     if (kk==2097)&&(ttt1(116,116)>0.01)
%         disp(kk);
%     end
        
%     if ( Energy(IDX_sorted(kk)) == 0 )
    % If the block energy < Th, then stop the matching pursuit. 
    if ( Energy(IDX_sorted(kk)) <= Th )
        break;
    else
        LocationMax = IDX_sorted(kk);
    end

    % get the block type and position that generates the maximum energy
    % ind2sub is slow. Do it manually.
%     [CoordinateMax,TypeBlk] = ind2sub([N,BNumber],LocationMax);
%     [CoordinateMax1,CoordinateMax2] = ind2sub([N1,N2],CoordinateMax);

    CoordinateMax = mod(LocationMax, N);
    if ( CoordinateMax == 0 )
        CoordinateMax = N;
    end
    TypeBlk = ceil(LocationMax / N);

    CoordinateMax1 = mod(CoordinateMax, N1);
    if ( CoordinateMax1 == 0 )
        CoordinateMax1 = N1;
    end
    CoordinateMax2 = ceil(CoordinateMax / N1);

    
    % set a flag on each the point in the selected block
    BlockChosen(CoordinateMax-1 + BGeo(1:BSize(TypeBlk),TypeBlk)) = 1;
    
    % Update of the energy by setting the energy of all the blocks that
    % intersect the selected block.
    % Only look in a local rectangle that includes the support the block
    % loop over block types.
    
    for k = 1:BNumber      
        % search in the local rectangle
        for i1 = max(CoordinateMax1-Support1,1) : min(CoordinateMax1+Support1,N1)     
            for i2 = max(CoordinateMax2-Support2,1) : min(CoordinateMax2+Support2,N2)
                n = i1+(i2-1)*N1;
                if (BPosition(n,k)==1)
                   % loop over all the points in the block being examined
                   for p = 1:BSize(k) 
                      q = n-1+BGeo(p,k);
                      % if it intersects with the selected block, then set
                      % its energy to 0
                      if (BlockChosen(q)==1)
%                           if n==29556
%                               disp('error');
%                           end
                            Energy(n,k) = 0;
                            break;
                      end
                   end                   
                end
            end
        end
    end
            
    BListSize = BListSize+1;
    BList(BListSize,:) = [CoordinateMax,TypeBlk];         
    Bepsilon(BListSize,:) = [CoordinateMax,Epsilon(LocationMax)];
    
    % reset the flag
    BlockChosen(CoordinateMax-1+BGeo(1:BSize(TypeBlk),TypeBlk)) = 0; 
           
    kk = kk + 1;
end

BList = BList(1:BListSize,:);
Bepsilon = Bepsilon(1:BListSize,:);
end