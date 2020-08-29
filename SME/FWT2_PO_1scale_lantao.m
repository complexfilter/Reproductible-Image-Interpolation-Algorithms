function wc = FWT2_PO_1scale_lantao(x,qmf)
% FWT2_PO -- 2-d MRA wavelet transform (periodized, orthogonal)
%  Usage
%    wc = FWT2_PO(x,L,qmf)
%  Inputs
%    x     2-d image (n by n array, n dyadic)
%    L     coarse level
%    qmf   quadrature mirror filter
%  Outputs
%    wc    2-d wavelet transform
%
%  Description
%    A two-dimensional Wavelet Transform is computed for the
%    array x.  To reconstruct, use IWT2_PO.
%
%  See Also
%    IWT2_PO, MakeONFilter
%
% [n,J] = quadlength(x);
[n1, n2] = size(x);
wc = x;
% nc = n;
% for jscal=J-1:-1:L,
top = (n2/2+1):n2; bot = 1:(n2/2);
for ix=1:n1
    row = wc(ix,1:n2);
    wc(ix,bot) = DownDyadLo(row,qmf);
    wc(ix,top) = DownDyadHi(row,qmf);
end
top = (n1/2+1):n1; bot = 1:(n1/2);
for iy=1:n2
    row = wc(1:n1,iy)';
    wc(top,iy) = DownDyadHi(row,qmf)';
    wc(bot,iy) = DownDyadLo(row,qmf)';
end
end
function d = DownDyadLo(x,qmf)
d = aconv(qmf,x);
n = length(d);
d = d(1:2:(n-1));
end
function d = DownDyadHi(x,qmf)
% DownDyadHi -- Hi-Pass Downsampling operator (periodized)
%  Usage
%    d = DownDyadHi(x,f)
%  Inputs
%    x    1-d signal at fine scale
%    f    filter
%  Outputs
%    y    1-d signal at coarse scale
%
%  See Also
%    DownDyadLo, UpDyadHi, UpDyadLo, FWT_PO, iconv
%
d = iconv( MirrorFilt(qmf),lshift(x));
n = length(d);
d = d(1:2:(n-1));
end

function y = MirrorFilt(x)
% MirrorFilt -- Apply (-1)^t modulation
%  Usage
%    h = MirrorFilt(l)
%  Inputs
%    l   1-d signal
%  Outputs
%    h   1-d signal with DC frequency content shifted
%        to Nyquist frequency
%
%  Description
%    h(t) = (-1)^(t-1)  * x(t),  1 <= t <= length(x)
%
%  See Also
%    DyadDownHi
%

y = -( (-1).^(1:length(x)) ).*x;
end


function y = iconv(f,x)
% iconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = iconv(f,x)
%  Inputs
%    f   filter
%    x   1-d signal
%  Outputs
%    y   filtered result
%
%  Description
%    Filtering by periodic convolution of x with f
%
%  See Also
%    aconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%
n = length(x);
p = length(f);
if p <= n
    xpadded = [x((n+1-p):n) x];
else
    z = zeros(1,p);
    for i=1:p
        imod = 1 + rem(p*n -p + i-1,n);
        z(i) = x(imod);
    end
    xpadded = [z x];
end
ypadded = filter(f,1,xpadded);
y = ypadded((p+1):(n+p));
end

function y = aconv(f,x)
% aconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = aconv(f,x)
%  Inputs
%    f    filter
%    x    1-d signal
%  Outputs
%    y    filtered result
%
%  Description
%    Filtering by periodic convolution of x with the
%    time-reverse of f.
%
%  See Also
%    iconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%

n = length(x);
p = length(f);
if p < n
    xpadded = [x x(1:p)];
else
    z = zeros(1,p);
    for i=1:p
        imod = 1 + rem(i-1,n);
        z(i) = x(imod);
    end
    xpadded = [x z];
end
fflip = f(end:-1:1);
ypadded = filter(fflip,1,xpadded);
y = ypadded(p:(n+p-1));
end
function y = lshift(x)
% lshift -- Circular left shift of 1-d signal
%  Usage
%    l = lshift(x)
%  Inputs
%    x   1-d signal
%  Outputs
%    l   1-d signal
%        l(i) = x(i+1) except l(n) = x(1)
%

y = [ x( 2:length(x) ) x(1) ];
end