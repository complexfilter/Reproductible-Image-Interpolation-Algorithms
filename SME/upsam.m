function y = upsam(x,S)
% y = upsam(x,S) inserts S-1 zeros between each term in the row vector x.
%  for example:   [1 0 2 0 3 0] = upsample([1 2 3]).   csb 3/1/93.
L = length(x);
y = [x;zeros(S-1,L)]; y = y(:).';
y = y(1:S*L-1);
end
