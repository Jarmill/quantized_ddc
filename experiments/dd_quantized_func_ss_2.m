rng(30, 'twister')

%a working example for superstability:

%more likely to be superstable when:
%   rho decreases to 1 (higher quantization density)
%   N_r increases (more buckets)
%   T increases (more data)

n = 3;
m = 2;

A0 = [-0.130044145839941,-0.397448693430848,0.202989496636112;-0.397448693430849,-0.499993279171062,0.299032700527649;0.202989496636112,0.299032700527648,-0.526192583398200];
B0 = [0.217864167978751,1.279977824637024;0.359211320717385,0;-1.155284985838844,0];

sys = struct('A', A0, 'B', B0);
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
T = 100;
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

B_r = 4;
% N_r = 11;
N_r = 9;
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
% rho = 1.3;
% rho = 1.5;
% rho = 1.4;
% rho = 1.05;
% rho = 0.7;


Nrho = 100;
rho_list = linspace(0.1, 1, Nrho);
lam_list = zeros(Nrho, 1);
% rho = 0.3;
for i = 1:length(rho_list)

% out = SS_quantized(sim, rho);
%     out = SS_quantized_sign(sim, rho_list(i));
    out = SS_quantized_clean_aff(rho_list(i), sys);
    if out.problem
        lam_list(i) = inf;
    else
        lam_list(i) = out.obj;
    end
end
% out = ESS_quantized_sign(sim, rho);

if out.problem==0
    K_rec = out.K;
    lam_rec = out.obj;
end


%%  plot the result
figure, 
hold on
plot(rho_list, lam_list_60', 'linewidth', 2)
plot(rho_list, lam_list_80', 'linewidth', 2)
plot(rho_list, lam_list_100', 'linewidth', 2)
plot(rho_list, lam_list_clean', 'linewidth', 2)
plot(xlim, [1, 1], '--k','linewidth', 2)
ylim([0, max(lam_list_80)])
legend({'T=60', 'T=80', 'T=100', 'truth', '\lambda=1'}, 'location', 'southwest')
xlabel('Quantization $\rho$', 'interpreter', 'latex')
ylabel('Gain $\lambda$', 'interpreter', 'latex')
title('Superstability Gain vs. Quantization', 'FontSize', 16)
