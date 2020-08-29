function imH = DirInterp2_3step(imL, m, n, itp_method)
%
% Directional image interpolation. Interpolate the image along the
% direcion arctan(m/n), with the interpolation kernel 'itp_method'. 
% 
% Input:
% imL: gray-level low-resolution image. It should be a 2D matrix. 
%
% m, n: theta=arctan(m/n) specifies the interpolation direction, with 
%   m>0, n>0 ==> 0 < theta < 90 degree
%   m>0, n<0 ==> 90 < theta < 180 degre
%   m = 0 ==> theta = 0 degre
%   n = 0 ==> theta = 90 degre
%
% itp_method specifies the interpolation kernel (see interp).
%   'spline': cubic spline.
%   'cubic': bicubic.
%   'linear': linear.
%
% Output:
% imH: directionally interpolated image. 
%
% % Reference:
% The directional interpolation is performed in 3 steps which are described
% in detail in the following paper. 
% S.Mallat and G.Yu, Super-Resolution with Sparse Mixing Estimators,
% submitted to IEEE Trans. on Image Processing, 2009.
% Website: http://www.cmap.polytechnique.fr/~mallat/SME
% 
% Contact: Guoshen Yu  yu@cmap.polytechnique.fr
%
% Copyright 2009 Guoshen Yu


[Nl1, Nl2] = size(imL);

% if the angle is horizontal ( m = 0 ) or vertical ( n = 0 ), do separable
% interpolation. (We can do nothing better anyway.)
if ( ( m == 0 ) || ( n == 0 ) )
    % linear interpolation
%     [CX,CY] = meshgrid(1:Nl1, 1:Nl2);
%     [X,Y] = meshgrid(1:0.5:Nl1, 1:0.5:Nl2);
    [CX,CY] = meshgrid(1:Nl2, 1:Nl1);
    [X,Y] = meshgrid(1:0.5:Nl2, 1:0.5:Nl1);
    imH = interp2(CX,CY,imL,X,Y, 'spline');
    return;
end

% if n < 0 which means that pi/2 < theta < pi, flip the image left-right
% and use -n instead of n. At the end of the algorith, reflip the
% high-resolution obtained.
flag_n_neg = 0;
if ( n < 0 )
    n = -n;
    flag_n_neg = 1;
    imL = fliplr(imL);
end

Nh1 = 2 * Nl1 - 1;
Nh2 = 2 * Nl2 - 1;

imH = zeros(Nh1, Nh2);

imH(1:2:end, 1:2:end) = imL;

% interprete the movement vector on the high-resolution image
mh = m * 2;
nh = n * 2;

% Step 1: interpolate line by line along theta
M_flag = zeros(Nh1, Nh2);
flag_finish = 0;
for i = 1 : 2 : Nh1
    for j = 1 : 2 : Nh2
        if ( M_flag(i,j) == 0 )      
            % trace the line
%             ind_ij = sub2ind([Nh1, Nh2], i, j);
            ind_ij = ( j - 1 ) * Nh1 + i;
            
            % the points along theta ABOVE (i,j)
            n1a = i-mh : -mh : 1;
            n2a = j+nh : nh : Nh2;
            l_min = min(length(n1a) ,length(n2a));
            n1a = n1a(1 : l_min);
            n2a = n2a(1 : l_min);

            % the points along theta BELOW (i,j)
            n1b = i+mh : mh : Nh1;
            n2b = j-nh : -nh : 1;
            l_min = min(length(n1b) ,length(n2b));
            n1b = n1b(1 : l_min);
            n2b = n2b(1 : l_min);
            
            if ( (~isempty(n1a)) && (~isempty(n2a)) )
%                 ind1 = sub2ind([Nh1, Nh2], n1a, n2a);
                ind1 = ( n2a - 1 ) * Nh1 + n1a;
            else
                ind1 = [];
            end

            % flip ind1 so that the line is in good order (from top to bottom)
            ind1 = fliplr(ind1);

            if ( (~isempty(n1b)) && (~isempty(n2b)) )
%                 ind2 = sub2ind([Nh1, Nh2], n1b, n2b);
                ind2 = ( n2b - 1 ) * Nh1 + n1b;
            else
                ind2 = [];
            end

            ind = [ind1 ind_ij ind2];

            % the crosses along theta
            line_LR = imH(ind);
    
            length_line = length(line_LR);
            
            % if only one cross on the line, skip
            if ( length_line == 1 ) 
                continue;
            end
            
            % Interpolation
            X = 1 : length_line;
            Xh = 1 : 0.5 : length_line;
            line_HR = interp1(X, line_LR, Xh, itp_method);
            
            % store the high-resolution line
%             [n1l, n2l] = ind2sub([Nh1, Nh2], ind);
            n1l = mod(ind, Nh1);
            n1l(n1l == 0) = Nh1;
            n2l = ceil(ind / Nh1);
            
            n1h = n1l(1) : m : n1l(end);
            n2h = n2l(1) : -n : n2l(end);
%             ind = sub2ind([Nh1, Nh2], n1h, n2h);
            ind = ( n2h - 1 ) * Nh1 + n1h;

            imH(ind) = line_HR;
            
            % set the flag
            M_flag(ind) = 1;
            
            % short cut if all have been visited
            tmp = M_flag(1:2:end, 1:2:end);
            if ( isempty(find(tmp==0, 1)) )
                flag_finish = 1;
                break;
            end
        end
    end
    if ( flag_finish == 1 )
        break;
    end
end
    

% Step 2: horizontal/vertical/diagonal interpolation
% m = 1, n even => vertical interpolation (because the vertical points have been filled in step 1)
% if ( ( m == 1 ) && ( mod(n,2) == 0 ) )
if ( ( mod(m,2) == 1 ) && ( mod(n,2) == 0 ) )
    % create an image n times big
    Nhh1 = (Nh1-1) * n + 1;
    Nhh2 = Nh2;
    imH2 = zeros(Nhh1, Nhh2);   
    imH2(1:n:end, :) = imH;
  
    X = 1 : n : Nhh1;
    Xh = 1 : Nhh1;
    for  j = 1 : 2 : Nh2
        % Interpolation
        line_LR = imH2(X, j);
        line_HR = interp1(X, line_LR, Xh, itp_method);
        imH2(Xh, j) = line_HR;
    end
    
% m even, n odd => horizontal interpolation (because the horizontal points have been filled in step 1)
% elseif ( ( mod(m,2) == 0) && ( n == 1 ) )
elseif ( ( mod(m,2) == 0) && ( mod(n,2) == 1 ) )
    % create an image n times big
    Nhh1 = Nh1;
    Nhh2 = (Nh2-1) * m + 1;
    imH2 = zeros(Nhh1, Nhh2);   
    imH2(:, 1:m:end) = imH;
     
    X = 1 : m : Nhh2;
    Xh = 1 : 1 : Nhh2;
    for  i = 1 : 2 : Nh1
        % Interpolation
        line_LR = imH2(i, X);
        line_HR = interp1(X, line_LR, Xh, itp_method);
        imH2(i, Xh) = line_HR;
    end
   
% m odd, n odd => diagonal interpolation (because the quincunx points have been filled in step 1)
% the required interpolation factor depends on m and n
% elseif ( (( m == 1 ) && ( n == 3 )) || (( m == 3 ) && ( n == 1 )) )
elseif ( ( ( m == 1 ) && ( mod(n,2) == 1 ) ) || ( ( mod(m,2) == 1 ) && ( n == 1 ) ) )
%     imH2 = zeros(Nh1*4-3, Nh2*4-3); 
%     [Nhh1, Nhh2] = size(imH2);
    % create an image n times big
    if ( m == 1 )
        delta = n + 1;
    else
        delta = m + 1;
    end
    Nhh1 = (Nh1-1) * delta + 1;
    Nhh2 = (Nh2-1) * delta + 1;
    imH2 = zeros(Nhh1, Nhh2);    
    imH2(1:delta:end, 1:delta:end) = imH; 
    % scan the diagonals (along 135 degree) and interpolation
    for  i = 1 : 2*delta : Nhh1        
        % trace the line
        n1l = i : delta : Nhh1;
        n2l = 1 : delta : Nhh2;        
        length_line = min(length(n1l) ,length(n2l));
        
        % if only one cross on the line, skip
        if ( length_line == 1 )
            continue;
        end
        
        n1l = n1l(1 : length_line);
        n2l = n2l(1 : length_line);
        
%         ind = sub2ind([Nhh1, Nhh2], n1l, n2l);
        ind = ( n2l - 1 ) * Nhh1 + n1l;
                        
        line_LR = imH2(ind);
                  
        % Interpolation
        X = 1 : length_line;
        Xh = 1 : 1/delta : length_line;

        line_HR = interp1(X, line_LR, Xh, itp_method);
        
        % store the high-resolution line       
        n1h = n1l(1) : n1l(end);
        n2h = n2l(1) : n2l(end);
%         ind = sub2ind([Nhh1, Nhh2], n1h, n2h);
        ind = ( n2h - 1 ) * Nhh1 + n1h;
        
        imH2(ind) = line_HR;
    end
    for  j = 2*delta+1 : 2*delta : Nhh2        
        % trace the line
        n1l = 1 : delta : Nhh1;
        n2l = j : delta : Nhh2;        
        length_line = min(length(n1l) ,length(n2l));
        
        % if only one cross on the line, skip
        if ( length_line == 1 )
            continue;
        end
        
        n1l = n1l(1 : length_line);
        n2l = n2l(1 : length_line);
        
%         ind = sub2ind([Nhh1, Nhh2], n1l, n2l);
        ind = ( n2l - 1 ) * Nhh1 + n1l;
                        
        line_LR = imH2(ind);
        
        % Interpolation
        X = 1 : length_line;
        Xh = 1 : 1/delta : length_line;

        line_HR = interp1(X, line_LR, Xh, itp_method);
        
        % store the high-resolution line       
        n1h = n1l(1) : n1l(end);
        n2h = n2l(1) : n2l(end);
%         ind = sub2ind([Nhh1, Nhh2], n1h, n2h);
        ind = ( n2h - 1 ) * Nhh1 + n1h;
        
        imH2(ind) = line_HR;
    end
else
    error('(m,n) is not defined in the dictionary!');
end


% Step 3: interpolate line by line along theta using the values obtained in
% step 2 to calculate the circles 
M_flag = zeros(Nhh1, Nhh2);
flag_finish = 0;
% m odd, n even 
if ( ( mod(m,2) == 1) && ( mod(n,2) == 0 ) )
% if ( ( m == 1 ) && ( mod(n,2) == 0 ) )
    % Scan over the points already interpolated
    % we only need to scan one column out of n. Others are systemmatically
    % covered.
    for i = 2 : 2 : Nhh1       
        for j = 1 : n : Nhh2   
            % trace the line
%             ind_ij = sub2ind([Nhh1, Nhh2], i, j);
            ind_ij = ( j - 1 ) * Nhh1 + i;
            
            % the points along theta ABOVE (i,j)
            n1a = i-2*m : -2*m : 1;
            n2a = j+2 : 2 : Nhh2;
            l_min = min(length(n1a) ,length(n2a));
            n1a = n1a(1 : l_min);
            n2a = n2a(1 : l_min);

            % the points along theta BELOW (i,j)
            n1b = i+2*m : 2*m : Nhh1;
            n2b = j-2 : -2 : 1;
            l_min = min(length(n1b) ,length(n2b));
            n1b = n1b(1 : l_min);
            n2b = n2b(1 : l_min);
            
            if ( (~isempty(n1a)) && (~isempty(n2a)) )
%                 ind1 = sub2ind([Nhh1, Nhh2], n1a, n2a);
                ind1 = ( n2a - 1 ) * Nhh1 + n1a;
            else
                ind1 = [];
            end

            % flip ind1 so that the line is in good order (from top to bottom)
            ind1 = fliplr(ind1);

            if ( (~isempty(n1b)) && (~isempty(n2b)) )
%                 ind2 = sub2ind([Nhh1, Nhh2], n1b, n2b);
                ind2 = ( n2b - 1 ) * Nhh1 + n1b;
            else
                ind2 = [];
            end

            ind = [ind1 ind_ij ind2];
            
             % if the line has been visited, skip
            if ( isempty(find(M_flag(ind)==0, 1))) 
                continue;
            end

            % the crosses along theta
            line_LR = imH2(ind);
    
            length_line = length(line_LR);
            
            % if only one cross on the line, skip
            if ( length_line == 1 ) 
                continue;
            end
            
            % Interpolation
            X = 1 : length_line;
            Xh = 1 : 0.5 : length_line;
            line_HR = interp1(X, line_LR, Xh, itp_method);
            
            % store the high-resolution line
%             [n1l, n2l] = ind2sub([Nhh1, Nhh2], ind);
            n1l = mod(ind, Nhh1);
            n1l(n1l == 0) = Nhh1;
            n2l = ceil(ind / Nhh1);
            
            n1h = n1l(1) : 1*m : n1l(end);
            n2h = n2l(1) : -1 : n2l(end);
%             ind = sub2ind([Nhh1, Nhh2], n1h, n2h);
            ind = ( n2h - 1 ) * Nhh1 + n1h;

            imH2(ind) = line_HR;
            
            % set the flag
            M_flag(ind) = 1;
            
            % short cut if all have been visited (except the boundary)
            tmp = M_flag(n+1:n:end-n, 2:2:end);
            if ( isempty(find(tmp==0, 1)) )
                flag_finish = 1;
                break;
            end            
        end
        
        if ( flag_finish == 1 )
            break;
        end
    end
       
    imH = imH2(1:n:end, :);  
% m even, n odd
% elseif ( ( mod(m,2) == 0) && ( n == 1 ) )
elseif ( ( mod(m,2) == 0) && ( mod(n,2) == 1 ) )
    % Scan over the points already interpolated
    % we only need to scan one column out of n. Others are systemmatically
    % covered.
    for i = 1 : m : Nhh1
        for j = 2 : 2 : Nhh2
            % trace the line
%             ind_ij = sub2ind([Nhh1, Nhh2], i, j);
            ind_ij = ( j - 1 ) * Nhh1 + i;
            
            % the points along theta ABOVE (i,j)
            n1a = i-2 : -2 : 1;
            n2a = j+2*n : 2*n : Nhh2;
            l_min = min(length(n1a) ,length(n2a));
            n1a = n1a(1 : l_min);
            n2a = n2a(1 : l_min);

            % the points along theta BELOW (i,j)
            n1b = i+2 : 2 : Nhh1;
            n2b = j-2*n : -2*n : 1;
            l_min = min(length(n1b) ,length(n2b));
            n1b = n1b(1 : l_min);
            n2b = n2b(1 : l_min);
            
            if ( (~isempty(n1a)) && (~isempty(n2a)) )
%                 ind1 = sub2ind([Nhh1, Nhh2], n1a, n2a);
                ind1 = ( n2a - 1 ) * Nhh1 + n1a;
            else
                ind1 = [];
            end

            % flip ind1 so that the line is in good order (from top to bottom)
            ind1 = fliplr(ind1);

            if ( (~isempty(n1b)) && (~isempty(n2b)) )
%                 ind2 = sub2ind([Nhh1, Nhh2], n1b, n2b);
                ind2 = ( n2b - 1 ) * Nhh1 + n1b;
            else
                ind2 = [];
            end

            ind = [ind1 ind_ij ind2];
            
            % if the line has been visited, skip
            if ( isempty(find(M_flag(ind)==0, 1))) 
                continue;
            end

            % the crosses along theta
            line_LR = imH2(ind);
    
            length_line = length(line_LR);
            
            % if only one cross on the line, skip
            if ( length_line == 1 ) 
                continue;
            end
            
            % Interpolation
            X = 1 : length_line;
            Xh = 1 : 0.5 : length_line;
            line_HR = interp1(X, line_LR, Xh, itp_method);
            
            % store the high-resolution line
%             [n1l, n2l] = ind2sub([Nhh1, Nhh2], ind);
            n1l = mod(ind, Nhh1);
            n1l(n1l == 0) = Nhh1;
            n2l = ceil(ind / Nhh1);
                        
            n1h = n1l(1) : 1 : n1l(end);
            n2h = n2l(1) : -1*n : n2l(end);
%             ind = sub2ind([Nhh1, Nhh2], n1h, n2h);
            ind = ( n2h - 1 ) * Nhh1 + n1h;

            imH2(ind) = line_HR;
            
            % set the flag
            M_flag(ind) = 1;
            
            % short cut if all have been visited 
            tmp = M_flag(2:2:end, n+1:n:end-n);
            if ( isempty(find(tmp==0, 1)) )
                flag_finish = 1;
                break;
            end            
        end
        
        if ( flag_finish == 1 )
            break;
        end
    end
       
    imH = imH2(:, 1:m:end);
% m odd, n odd
% elseif ( (( m == 1 ) && ( n == 3 )) || (( m == 3 ) && ( n == 1 )) )
elseif ( ( m == 1 ) && ( mod(n,2) == 1 ) )
    delta = n + 1;
    % scan only the values obtained in step 2
    % only points at one position (2+8*k1,2+8*k2) needs to be scanned since the
    % trace pass the points on the other positions (4+8*k1,4*8*k2), (6+8*k1,6*8*k2)
    % and (8+8*k1,8*8*k2)
    for i = 2 : 2*delta : Nhh1
        for j = 2 : 2*delta : Nhh2
            % trace the line
%             ind_ij = sub2ind([Nhh1, Nhh2], i, j);
            ind_ij = ( j - 1 ) * Nhh1 + i;
            
            % the points along theta ABOVE (i,j)
            n1a = i-2 : -2 : 1;
            n2a = j+2*n : 2*n : Nhh2;
            l_min = min(length(n1a) ,length(n2a));
            n1a = n1a(1 : l_min);
            n2a = n2a(1 : l_min);

            % the points along theta BELOW (i,j)
            n1b = i+2 : 2 : Nhh1;
            n2b = j-2*n : -2*n : 1;
            l_min = min(length(n1b) ,length(n2b));
            n1b = n1b(1 : l_min);
            n2b = n2b(1 : l_min);
            
            if ( (~isempty(n1a)) && (~isempty(n2a)) )
%                 ind1 = sub2ind([Nhh1, Nhh2], n1a, n2a);
                ind1 = ( n2a - 1 ) * Nhh1 + n1a;
            else
                ind1 = [];
            end

            % flip ind1 so that the line is in good order (from top to bottom)
            ind1 = fliplr(ind1);

            if ( (~isempty(n1b)) && (~isempty(n2b)) )
%                 ind2 = sub2ind([Nhh1, Nhh2], n1b, n2b);
                ind2 = ( n2b - 1 ) * Nhh1 + n1b;
            else
                ind2 = [];
            end

            ind = [ind1 ind_ij ind2];
            
            % if the line has been visited, skip
            if ( isempty(find(M_flag(ind)==0, 1)))
                continue;
            end

            % the crosses along theta
            line_LR = imH2(ind);
    
            length_line = length(line_LR);
            
            % if only one cross on the line, skip
            if ( length_line == 1 ) 
                continue;
            end
            
            % Interpolation
            X = 1 : length_line;
            Xh = 1 : 0.5 : length_line;
            line_HR = interp1(X, line_LR, Xh, itp_method);
            
            % store the high-resolution line
%             [n1l, n2l] = ind2sub([Nhh1, Nhh2], ind);
            n1l = mod(ind, Nhh1);
            n1l(n1l == 0) = Nhh1;
            n2l = ceil(ind / Nhh1);
            
            n1h = n1l(1) : 1 : n1l(end);
            n2h = n2l(1) : -n : n2l(end);
%             ind = sub2ind([Nhh1, Nhh2], n1h, n2h);
            ind = ( n2h - 1 ) * Nhh1 + n1h;

            imH2(ind) = line_HR;
            
            % set the flag
            M_flag(ind) = 1;                                
        end        
    end
       
    imH = imH2(1:delta:end, 1:delta:end);
    
elseif ( ( mod(m,2) == 1 ) && ( n == 1 ) )
    delta = m + 1;
    % scan only the values obtained in step 2
    % only points at one position (2+8*k1,2+8*k2) needs to be scanned since the
    % trace pass the points on the other positions (4+8*k1,4*8*k2), (6+8*k1,6*8*k2)
    % and (8+8*k1,8*8*k2)
    for i = 2 : 2*delta : Nhh1
        for j = 2 : 2*delta : Nhh2
            % trace the line
%             ind_ij = sub2ind([Nhh1, Nhh2], i, j);
            ind_ij = ( j - 1 ) * Nhh1 + i;
            
            % the points along theta ABOVE (i,j)
            n1a = i-2*m : -2*m : 1;
            n2a = j+2 : 2 : Nhh2;
            l_min = min(length(n1a) ,length(n2a));
            n1a = n1a(1 : l_min);
            n2a = n2a(1 : l_min);

            % the points along theta BELOW (i,j)
            n1b = i+2*m : 2*m : Nhh1;
            n2b = j-2 : -2 : 1;
            l_min = min(length(n1b) ,length(n2b));
            n1b = n1b(1 : l_min);
            n2b = n2b(1 : l_min);
            
            if ( (~isempty(n1a)) && (~isempty(n2a)) )
%                 ind1 = sub2ind([Nhh1, Nhh2], n1a, n2a);
                ind1 = ( n2a - 1 ) * Nhh1 + n1a;
            else
                ind1 = [];
            end

            % flip ind1 so that the line is in good order (from top to bottom)
            ind1 = fliplr(ind1);

            if ( (~isempty(n1b)) && (~isempty(n2b)) )
%                 ind2 = sub2ind([Nhh1, Nhh2], n1b, n2b);
                ind2 = ( n2b - 1 ) * Nhh1 + n1b;
            else
                ind2 = [];
            end

            ind = [ind1 ind_ij ind2];
            
            % if the line has been visited, skip
            if ( isempty(find(M_flag(ind)==0, 1)))
                continue;
            end

            % the crosses along theta
            line_LR = imH2(ind);
    
            length_line = length(line_LR);
            
            % if only one cross on the line, skip
            if ( length_line == 1 ) 
                continue;
            end
            
            % Interpolation
            X = 1 : length_line;
            Xh = 1 : 0.5 : length_line;
            line_HR = interp1(X, line_LR, Xh, itp_method);
            
            % store the high-resolution line
%             [n1l, n2l] = ind2sub([Nhh1, Nhh2], ind);
            n1l = mod(ind, Nhh1);
            n1l(n1l == 0) = Nhh1;
            n2l = ceil(ind / Nhh1);
            
            n1h = n1l(1) : m : n1l(end);
            n2h = n2l(1) : -1 : n2l(end);
%             ind = sub2ind([Nhh1, Nhh2], n1h, n2h);
            ind = ( n2h - 1 ) * Nhh1 + n1h;

            imH2(ind) = line_HR;
            
            % set the flag
            M_flag(ind) = 1;                                
        end        
    end
       
    imH = imH2(1:delta:end, 1:delta:end);    
else
    error('(m,n) is not defined in the dictionary!');
end
    
if ( flag_n_neg == 1 )
    imH = fliplr(imH);
end
