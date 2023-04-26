rng(20, 'twister')
n = 3;
m = 2;
% ss = rss(n, n, m);
A0 = 0.9*eye(3)+ 0.2*ones(3);
% A0(3, 2) = 2;
B0 = [1 0; 0 1; 0 1];
ss = struct('A', A0, 'B', B0)
% ss.A = 1.5*ss.A;

%perform logarithmic quantization with level q
% q = 0.1;
q = 0.1665;
rho = (1+q)/(1-q);

v = sdpvar(n, 1);
S = sdpvar(m, n);

Y = diag(v);

% Acl = ss.A+ss.B*K;
Acl_low = ss.A*Y + ss.B*S*(1-q);
Acl_high = ss.A*Y + ss.B*S*(1+q);
Acl = ss.A*Y + ss.B*S;
M = sdpvar(n, n);
delta = 1e-2;
cons = [(M(:)- Acl_low(:)) >=0 ; (M(:)+Acl_low(:))>= 0;...
    (M(:)- Acl_high(:)) >=0 ; (M(:)+Acl_high(:))>= 0;...
    sum(M, 2) <= v-delta; v>=delta];
opts = sdpsettings;
sol = optimize(cons, [], opts);

disp(sol.info)

%% recovery
Mr = value(M);
Sr = value(S);
vr = value(v);
Aclr = value(Acl);
Kr = Sr*diag(1./vr);
Aclr_c = Aclr*diag(1./vr);


