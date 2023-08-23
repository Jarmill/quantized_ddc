
rng(60, 'twister')
% n = 4;
% m = 3;

n = 3;
m = 2;
ss = drss(n, n, m);
% ss.A = ss.A + 0.3*eye(n);
% eig(ss.A)
ss.A = ss.A*2;
% ss.A = ss.A*1.4;
% ss.A = ss.A*1.3;
% rng(20, 'twister')
% n = 3;
% m = 2;
% % ss = rss(n, n, m);
% A0 = 0.9*eye(3)+ 0.2*ones(3);
% % A0(3, 2) = 2;
% B0 = [1 0; 0 1; 0 1];
% ss = struct('A', A0, 'B', B0)
% ss.A = 1.5*ss.A;

%perform logarithmic quantization with level q
% q = 0.1;
rho = 0.45;
% rho = 0.9;
q = (1-rho)/(1+rho);

v = sdpvar(n, 1);
S = sdpvar(m, n);

Y = diag(v);

% Acl = ss.A+ss.B*K;
Acl_11= ss.A*Y + ss.B*(eye(2)+diag([1;1])*q)*S;
Acl_12= ss.A*Y + ss.B*(eye(2)+diag([1;-1])*q)*S;
Acl_21= ss.A*Y + ss.B*(eye(2)+diag([-1;1])*q)*S;
Acl_22= ss.A*Y + ss.B*(eye(2)+diag([-1;-1])*q)*S;
Acl = ss.A*Y + ss.B*S;
M = sdpvar(n, n);
delta = 1e-2;
cons = [(M(:)- Acl_11(:)) >=0 ; (M(:)+Acl_11(:))>= 0;...
    (M(:)- Acl_12(:)) >=0 ; (M(:)+Acl_12(:))>= 0;...
    (M(:)- Acl_21(:)) >=0 ; (M(:)+Acl_21(:))>= 0;...
    (M(:)- Acl_22(:)) >=0 ; (M(:)+Acl_22(:))>= 0;...    
    sum(M, 2) <= v-delta; v>=delta];
opts = sdpsettings;
sol = optimize(cons, [], opts);

disp(sol.info)

%% recovery
Mr = value(M);
Sr = value(S);
vr = value(v);
Aclv = value(Acl); %v-weighted closed loop system
Kr = Sr*diag(1./vr);
Aclr = Aclv*diag(1./vr); %nominal closed-loop system


%% simulate a trajectory

