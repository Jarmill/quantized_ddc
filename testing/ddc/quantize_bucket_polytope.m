rng(50, 'twister')
n = 4;
m = 2;
[imat, jmat] = meshgrid(1:n, 1:n);

A0 = (1/n)*min(imat./jmat, jmat./imat) + 0.5*eye(n);
B0 = eye(n, m);


%% sample data

% T = 60;
% T = 80;
% T = 100;
% T = 120;
T = 200;
% T = 350;
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

%% this code will go into data_cons_bucket

%get bounds from the buckets
lb = zeros(n, T);
ub = zeros(n, T);
for i = 1:Nbucket
    lb(Sb{i}) = buckets(i, 1);
    ub(Sb{i}) = buckets(i, 2);
end

mask_lb_inf = (lb~=-Inf);
mask_ub_inf = (ub~=Inf);
%form expressions
CA0 = kron(Xn', speye(n));
CB0 = kron(U', speye(n));
d0 = reshape(Xp, [], 1);

%extract the non-infinite bounds
CA_ub = CA0(mask_ub_inf(:), :);
CB_ub = CB0(mask_ub_inf(:), :);

CA_lb = CA0(mask_lb_inf(:), :);
CB_lb = CB0(mask_lb_inf(:), :);

% CA = [CA_ub; -CA_lb];
% CB = [CB_ub; -CB_lb];

C = [CA_ub, CB_ub; -CA_lb, -CB_lb];
d = [ub(mask_ub_inf(:)); -lb(mask_lb_inf(:))];

%reduce the dimension
[C_red, d_red] = nontrivial_constraints(C, d);

%% check validity
sys_true = [A0(:); B0(:)];
res_true = d - C*sys_true;
res_red_true = d_red - C_red*sys_true;