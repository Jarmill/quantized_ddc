rng(20, 'twister')
n = 3;
m = 2;
% ss = rss(n, n, m);
% A0 = 0.6*eye(3)- 0.1*ones(3);
A0 = 0.9*eye(3)+ 0.2*ones(3);
% A0(3, 2) = 2;
B0 = [1 0; 0 1; 0 1];
ss = struct('A', A0, 'B', B0);

%% sample data
% T = 120;
% T = 300;
T = 400;
% T = 90;
umax = 1;
Xn = 2*rand(n, T)-1;
U = umax*(2*rand(m, T)-1);
Xp = A0*Xn + B0*U;


buckets = [-inf, -1; -1, 0; 0, 1; 1, inf];
Nbucket = size(buckets, 1);
Sb = cell(Nbucket, 1);
for i = 1:Nbucket
    Sb{i} = (Xp >= buckets(i, 1))  & (Xp <= buckets(i, 2));
end



%% synthesize a controller
% ss.A = 1.5*ss.A;

%perform logarithmic quantization with level q
% q = 0.1;
% q = 0.1;
% q = 0.05;
q = 0.04;
rho = (1+q)/(1-q);

% K = sdpvar(m, n, 'full');
A = sdpvar(n, n, 'full');
B = sdpvar(n, m, 'full');

v = sdpvar(n, 1);
S = sdpvar(m, n);

Y = diag(v);

Xp_curr = A*Xn + B*U;
Xp_sub = replace(Xp_curr, [A(:); B(:)], [A0(:); B0(:)]);

rob_vars = [A(:); B(:)];
M0 = sdpvar(n, n);
Ma = cell(length(rob_vars), 1);

%affine adjustable robust counterpart
M = M0;
for i = 1:length(rob_vars)
    Ma{i} = sdpvar(n, n);
    M = M +  Ma{i}*rob_vars(i);
end


dd_term = [];
for i = 1:Nbucket
    if buckets(i, 1) > -inf
        dd_low = Xp_curr - buckets(i, 1);
        dd_term = [dd_term; dd_low(Sb{i})];
    end
    if buckets(i, 2) < inf
        dd_high = (buckets(i, 2) - Xp_curr);
        dd_term = [dd_term; dd_high(Sb{i})];
    end
end
dd_con = (dd_term>=0);

%(extended super) stabilization
% Acl = A+B*K;
% Acl = A*Y + B*S;
Acl_1  = A*Y + B*(eye(2) - diag([1, 1])*q)*S;
Acl_2  = A*Y + B*(eye(2) - diag([1, -1])*q)*S;
Acl_3  = A*Y + B*(eye(2) - diag([-1, 1])*q)*S;
Acl_4  = A*Y + B*(eye(2) - diag([-1, -1])*q)*S;

delta = 1e-2;

ss_con = [(M(:)-Acl_1(:)) >=0 ; (M(:)+Acl_1(:))>= 0; ...
    (M(:)-Acl_2(:)) >=0 ; (M(:)+Acl_2(:))>= 0; ...
    (M(:)-Acl_3(:)) >=0 ; (M(:)+Acl_3(:))>= 0; ...
    (M(:)-Acl_4(:)) >=0 ; (M(:)+Acl_4(:))>= 0; ...    
    sum(M, 2) <= (v-delta)];
v_con = [v>=delta; sum(v)==1];

% M = sdpvar(n, n);

cons = [v_con:'lyapunov'; ss_con:'superstability'; dd_con:'buckets'; uncertain(rob_vars)];
opts = sdpsettings('robust.lplp', 'duality');
sol = optimize(cons, [], opts);

disp(sol.info)

%% recovery

M_rep = replace(M, rob_vars, [A0(:); B0(:)]);
Mr = value(M_rep);
Sr = value(S);
vr = value(v);
% Aclr = value(Acl);
Aclr0 = value(A0 + B0*K);
Aclr_low = value(A0 + (1-delta)*B0*Kr);
Aclr_high = value(A0 + (1+delta)*B0*Kr);
