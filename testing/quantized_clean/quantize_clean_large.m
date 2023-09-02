rng(20, 'twister')

%rho = 0.5;
n = 8;
m = 4;

%rho = 0.6;
n = 5;
m = 3;



% n = 14;
% m = 7;
% ss = rss(n, n, m);
% A0 = eye(n)*1-ones(n)*0.2;
[imat, jmat] = meshgrid(1:n, 1:n);

A0 = (1/n)*min(imat./jmat, jmat./imat) + 0.45*eye(n);

% A0 = -0.8*eye(n) + ones(n)*0.1;
% A0 = 1./gallery('minij',n);
% A0(3, 2) = 2;
B0 = eye(n, m);
% B0 = B0(end:-1:1,end:-1:1 );
% e = ones(m, 1);
% B0 = full(spdiags([-e e e],-1:1,n,m));
% A0 = eye(3);
% B0 = [0 0 ; 1 0; 0 1];
ss = struct('A', A0, 'B', B0);

rho = 0.25;
% rho = 0.001;

% out = ESS_quantized_sign([], rho, ss);
% out = ESS_quantized_clean_aff(rho, ss);
out = ESS_quantized_clean_aff(rho, ss);


%% recovery
if ~out.problem

    Acl_out = ss.A + ss.B * out.K;
    [eig(Acl_out), eig(A0)]'
else
    eig(A0)'
end
%% simulate a trajectory
 