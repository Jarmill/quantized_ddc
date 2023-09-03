rng(30, 'twister')

%a working example for superstability:

%more likely to be superstable when:
%   rho decreases to 1 (higher quantization density)
%   N_r increases (more buckets)
%   T increases (more data)

%large model runs out of memory in MATLAB
%requires on the order of 432 GB
% n = 14;
% m = 7;


%medium size
% n = 8;
% m = 4;

n=5;
m=3;
[imat, jmat] = meshgrid(1:n, 1:n);

A0 = (1/n)*min(imat./jmat, jmat./imat) + 0.45*eye(n);
B0 = eye(n, m);

% ss = drss(n, n, m);
% A0 = ss.A + 0.1*eye(n);
% ss.A = ss.A*1.4;
% A0 = 0.6*eye(3)- 0.1*ones(3);

% A0 = ss.A;
% B0 = ss.B;

% A0 = 0.9*eye(3)+ 0.2*ones(3);
% A0 = diag([1.2; 0.5; 0.5]) - 0.1*ones(3);
% A0(3, 2) = 2;
% B0 = [1 0; 0 1; 0 1];
% B0 = [0 1; 0 1; 1 0];
% ss = struct('A', A0, 'B', B0);


%% sample data

% T = 60;
% T = 80;
% T = 100;
% T = 120;
% T = 200;
% T = 300;
T = 350;
% T = 400;
% T = 450;
% T = 90;
umax = 10;
Xn = 2*rand(n, T)-1;
U = umax*(2*rand(m, T)-1);
Xp = A0*Xn + B0*U;


% buckets = [-inf, -1; -1, 0; 0, 1; 1, inf];
% buckets = 


%n=5, m=3, T = 350, rho=1. works
% rho = 1;
B_r = 6;
N_r = 25;


% B_r = 8;
% N_r = 17;
% N_r = 21;
% B_r = 4;
% N_r = 11;
% N_r = 9;
% N_r = 6;


% B_r = 6;
% N_r = 9;

fin_range = linspace(-B_r, B_r, N_r);
buckets = [[-inf, fin_range]; [fin_range, inf]]'; 
Nbucket = size(buckets, 1);
Sb = cell(Nbucket, 1);
for i = 1:Nbucket
    Sb{i} = (Xp >= buckets(i, 1))  & (Xp <= buckets(i, 2));
end


sim = struct('X', Xn, 'U', U, 'buckets', buckets);
sim.Sb =  Sb;

%run the quantizer

%rho should be less than 1?

% rho = 0.75;
% rho = 1;
rho = 0.8;

% out = ESS_quantized(sim, rho);
out = ESS_quantized_sign(sim, rho);


% Nrho = 100;
% rho_list = linspace(0.1, 1, Nrho);
% lam_list = zeros(Nrho, 1);
% % rho = 0.3;
% for i = 1:length(rho_list)
% 
% % out = SS_quantized(sim, rho);
%     out = SS_quantized_sign(sim, rho_list(i));
%     if out.problem
%         lam_list(i) = inf;
%     else
%         lam_list(i) = out.obj;
%     end
% end
% % out = ESS_quantized_sign(sim, rho);
% 
% if out.problem==0
%     K_rec = out.K;
%     lam_rec = out.obj;
% end
