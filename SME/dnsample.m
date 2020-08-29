function y = dnsample(x)
% y = dnsample(x) samples x by removing the even terms in x.
%  for example:   [1 3] = dnsample([1 2 3 4]).   csb 3/1/93.
L = length(x);
y = x(1:2:L);
end