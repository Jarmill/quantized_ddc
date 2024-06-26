function [out] = SS_quantized(sim, rho)
%SS_Quantized superstabilizing control for quantized linear systems
% use affine-adjustable robust counterparts

%process the input

n = size(sim.X, 1);
m = size(sim.U, 1);

if length(rho) == 1
    rho = ones(m, 1)*rho;
end

q = (1-rho)./(1+rho);


%declare the variables

%uncertain (temporary)
A = sdpvar(n, n, 'full');
B = sdpvar(n, m, 'full');
rob_vars = [A(:); B(:)];

%design parameters
K = sdpvar(m, n, 'full');
lambda = sdpvar(1, 1);

M0 = sdpvar(n, n);
Ma = cell(length(rob_vars), 1);

%affine adjustable robust counterpart
M = M0;
for i = 1:length(rob_vars)
    Ma{i} = sdpvar(n, n);
    M = M +  Ma{i}*rob_vars(i);
end

signs = 2*ff2n(m)-1;


%start up constraint generation
lam_cons = [sum(M, 2) <= lambda];

%iterate over all sign patterns
I = eye(m);
ss_cons = [];
for i = 1:size(signs, 1)
    sign_curr = signs(i, :);
    
    Delta = diag(sign_curr'.*q);
    
    Acl_curr = A + B*(Delta + I)*K;
    
    ss_cons = [ss_cons; (M(:)-Acl_curr (:)) >=0 ; (M(:)+Acl_curr (:))>= 0];
end

%data processing

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


cons = [dd_con:'data'; lam_cons:'objective'; ...
    ss_cons:'M bounding'; uncertain(rob_vars)];

%solve the program
opts = sdpsettings('robust.lplp', 'duality');
sol = optimize(cons, lambda, opts);

disp(sol.info)

out = struct;
out.problem = sol.problem;
if sol.problem ==0
    out.K = value(K);
    out.obj = value(lambda);
    out.lambda = out.obj;
    out.M0 = value(M0);
    out.Ma = cellfun(@(m) value(m), Ma, 'UniformOutput', false);
    out.info = sol.info;    
    out.problem = sol.problem;
    disp(out.obj)
end


end

