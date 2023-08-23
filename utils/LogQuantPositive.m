% This routine is from the example in the online tutorial:
% https://www.mathworks.com/help/control/ug/absolute-stability-for-quantized-system.html
% with supporting code 

% See the reference paper "The sector bound approach to quantized feedback
% control," Minyue Fu and Lihua Xie, IEEE Transactions on Automatic Control
% 50(11), 2005, 1698-1711.

function y = LogQuantPositive(u,u0,rho)


% Quantization level belongs to (alpha,beta] and beta-alpha = 1. Note that
% beta = log((1-delta)*u/u0)/log(rho);
delta = (1-rho)/(1+rho);
alpha = log((1+delta)*u/u0)/log(rho);

% Level 
level = ceil(alpha);

% Value after quantization
y = rho^level*u0;


