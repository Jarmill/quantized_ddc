rng(20, 'twister')
n = 3;
m = 2;
% ss = rss(n, n, m);
% A0 = 0.6*eye(3)- 0.1*ones(3);

% A0 = 0.9*eye(3)+ 0.2*ones(3);
A0 = diag([1.2; 0.5; 0.5]) - 0.1*ones(3);
% A0(3, 2) = 2;
B0 = [1 0; 0 1; 0 1];
% B0 = [0 1; 0 1; 1 0];
ss = struct('A', A0, 'B', B0);

%% sample data

% T = 60;
T = 80;
% T = 120;
% T = 200;
% T = 400;
% T = 90;
umax = 10;
Xn = 2*rand(n, T)-1;
U = umax*(2*rand(m, T)-1);
Xp = A0*Xn + B0*U;


% buckets = [-inf, -1; -1, 0; 0, 1; 1, inf];
% buckets = 

% B_r = 4;
% N_r = 11;
% N_r = 6;


B_r = 6;
N_r = 9;

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

rho = 0.8;

% out = SS_quantized(sim, rho);
% out = SS_quantized_sign(sim, rho);
out = ESS_quantized_sign(sim, rho);

if out.problem==0
    K_rec = out.K;
    lam_rec = out.obj;
end