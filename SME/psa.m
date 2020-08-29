function p = psa(h,kk)
% p = psa(h,kk)   calculates samples of the scaling function
%  phi(t) = p  by  kk  successive approximations from the
%  scaling coefficients  h.  Initial iteration is a constant.
%  phi_k(t) is plotted at each iteration.     csb 5/19/93
%
if nargin==1
    kk=11; 
end                              % Default number of iterations
h2= h*2/sum(h);                  % normalize  h(n)
K = length(h2)-1; S = 1280;       % Sets sample density
p = [ones(1,3*S*K),0]/(3*K);     % Sets initial iteration
P = p(1:K*S);                    % Store for later plotting
%axis([0 K*S+2 -.5 1.4]);
hu = upsam(h2,S);                % upsample h(n) by S
for iter = 0:kk                  % Successive approx.
   p = dnsample(conv(hu,p));     % convolve and down-sample
   %plot(p); pause;               % plot each iteration
%   P = [P;p(1:K*S)];            % store each iter. for plotting
end
p = p(1:K*S);                    % only the supported part
% L = length(p);
% x = ([1:L])/(S);
% axis([0 3 -.5 1.4]);
% plot(x,p);                       % Final plot
% title('Scaling Function by Successive Approx.');
% ylabel('Scaling Function');
% xlabel('x');

end