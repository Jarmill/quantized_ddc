rng(20, 'twister')
n = 3;
m = 2;
% ss = rss(n, n, m);
A0 = 0.8*eye(3)+ 0.2*ones(3);
% A0(3, 2) = 2;
B0 = [1 0; 0 1; 0 1];
% A0 = eye(3);
% B0 = [0 0 ; 1 0; 0 1];
ss = struct('A', A0, 'B', B0)
% ss.A = 1.5*ss.A;

K = sdpvar(m, n);

Acl = ss.A+ss.B*K;
M = sdpvar(n, n);
delta = 1e-2;
cons = [M(:)- Acl(:) >=0 ; M(:) + Acl(:)>= 0; sum(M, 2) <= (1-delta)];
opts = sdpsettings;
sol = optimize(cons, [], opts);

%% recovery
Mr = value(M);
Kr = value(K);
Aclr = value(Acl);

%% simulate a trajectory
