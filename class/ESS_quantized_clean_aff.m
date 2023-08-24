function [out] = ESS_quantized_clean_aff(rho, sys)
%ESS_Quantized extended superstabilizing control for quantized linear systems
% use affine-adjustable robust counterparts

%process the input

[n, m] = size(sys.B);

A = sys.A;
B = sys.B;
if length(rho) == 1
    rho = ones(m, 1)*rho;
end

q = (1-rho)./(1+rho);


%declare the variables

%design parameters
v = sdpvar(n, 1);
S = sdpvar(m, n);

Y = diag(v);
lambda = sdpvar(1, 1);

M0 = sdpvar(n, n);

%affine adjustable robust counterpart
M = M0;

signs = 2*ff2n(m)-1;


%start up constraint generation
delta = 1e-3;
lam_cons = [sum(M, 2) <= (v - delta)];

%iterate over all sign patterns
I = eye(m);
ss_cons = [];
for i = 1:size(signs, 1)
    sign_curr = signs(i, :);
    
    Delta = diag(sign_curr'.*q);
    
    Acl_curr = A*Y + B*(Delta + I)*S;
    
    ss_cons = [ss_cons; (M(:)-Acl_curr (:)) >=0 ; (M(:)+Acl_curr (:))>= 0];
end

%data processing


v_con = [v>=delta; lambda >= v];
cons = [v_con:'v';  lam_cons:'objective'; ...
    ss_cons:'M bounding'];

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
    out.obj = 0;
    
    out.M0 = value(M0);
%     out.Ma = cellfun(@(m) value(m), Ma, 'UniformOutput', false);
    out.info = sol.info;       
    disp(out.obj)    
end


end

