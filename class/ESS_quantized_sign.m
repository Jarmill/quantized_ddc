function [out] = ESS_quantized_sign(sim, rho, sys)
%SS_Quantized superstabilizing control for quantized linear systems
% use affine-adjustable robust counterparts

%process the input

if nargin < 3
    n = size(sim.X, 1);
    m = size(sim.U, 1);
else
    [n, m] = size(sys.B); 
end

if length(rho) == 1
    rho = ones(m, 1)*rho;
end

q = (1-rho)./(1+rho);


%declare the variables

%uncertain (temporary)
if nargin < 3
    ROBUST = 1;
    A = sdpvar(n, n, 'full');   
    B = sdpvar(n, m, 'full');
    rob_vars = [A(:); B(:)];
else
    ROBUST = 0;
    A = sys.A;
    B = sys.B;
end


%design parameters
v = sdpvar(n, 1);
S = sdpvar(m, n);
Y = diag(v);

lambda = sdpvar(1, 1);


signs_u = 2*ff2n(m)-1;

signs_x = 2*ff2n(n) - 1;


%iterate over all sign patterns
I = eye(m);
ss_cons = [];
for i = 1:size(signs_u, 1)
    sign_curr = signs_u(i, :);
    
    Delta = diag(sign_curr'.*q);
    
    Acl_curr = A*Y + B*(Delta + I)*S;
    
    for j = 1:size(signs_x, 1)
        sign_x_curr = signs_x(j, :);
        Acl_sign = Acl_curr*sign_x_curr';
        ss_cons = [ss_cons; v >= Acl_sign ];
    end
    
end

%data processing

if ROBUST
Xp_curr = A*sim.X + B*sim.U;

dd_term = [];
Nbucket = size(sim.buckets, 1);
for i = 1:Nbucket
    if sim.buckets(i, 1) > -inf
        dd_low = Xp_curr - sim.buckets(i, 1);
        dd_term = [dd_term; dd_low(sim.Sb{i})];
    end
    if sim.buckets(i, 2) < inf
        dd_high = (sim.buckets(i, 2) - Xp_curr);
        dd_term = [dd_term; dd_high(sim.Sb{i})];
    end
end
dd_con = (dd_term>=0);

rob_con = uncertain(rob_vars);
else
    dd_con = [];
    rob_con = [];
end

cons = [dd_con; lambda>=v; v>=1e-3;sum(v)==1;...
    ss_cons:'Sign-Superstability'; rob_con];

%solve the program
opts = sdpsettings('robust.lplp', 'duality');
sol = optimize(cons, lambda, opts);

disp(sol.info)

out = struct;
out.problem = sol.problem;
if sol.problem ==0
    out.S = value(S);
    out.v = value(v);
    out.K = out.S*diag(1./out.v);
    out.obj = value(lambda);
    out.lambda = out.obj;
    out.info = sol.info;    
    out.problem = sol.problem;
    disp(out.obj)
end


end

