function [wi, fxi] = RMLS_mm(xi, theta, Y_lab, vec_G_lab, lamda, S, Y_hat,  deta, w_hat, gama)
%RMLS Summary of this function goes here
%   Detailed explanation goes here
% This function is the solution. 
[dim, num_lab] = size( vec_G_lab);
D_theta = diag(theta);
J1 = Y_lab* D_theta * vec_G_lab';

lamda = trace(vec_G_lab * D_theta * vec_G_lab')/10^6;

M1 = vec_G_lab * D_theta * vec_G_lab'+ lamda * eye(dim,dim);



if(nargin < 6)
    wi = J1 * inv(M1);
    fxi = wi * xi;
else
      D_S = diag(S);
      [ddim, num_i] = size(Y_hat);
      X = repmat(xi,[1,num_i]);
      J2 = Y_hat * D_S * X';
      M2= (S'*ones(num_i,1))*(xi*xi');

%     [r,w,n] = size(w_hat);
%     ww_hat = reshape(w_hat, r*w,n);
%     ww_hat = ww_hat* S;
%     J3 = reshape(ww_hat,r,w);
%     M3 = ( S'*ones(num_i,1)) * eye(dim,dim);
%      wi = (J1 + deta*J2 + gama*J3)* inv(M1 + deta*M2 + gama*M3);
      wi = (J1 + deta*J2 )* inv(M1 + deta*M2);
%     wi = (J1  + gama*J3)* inv(M1  + gama*M3);
%      wi = J1 * inv(M1);
      fxi = wi * xi;
end
end
