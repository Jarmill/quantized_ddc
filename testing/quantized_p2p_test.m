
%try out the peak-to-peak formulation of quantized DDC on a clean system
%if the objective is less than 1, then we should have stabilization over
%the quantized system.
%
%need to check the exact gain formulas though
rng(20, 'twister')
n = 3;
m = 2;
ss = rss(n, n, m);

K = sdpvar(m, n);


lam_ref = 0.6075;
lam_test = 0.7;

% delta = 

%feasible: lam_test = 0.7, rho = 0.98
% rho = 0.98;
% rho = 0.95;
rho = 0.1;
delta = (1+rho)/(1-rho);

nu = norm(ss.B/delta, 'inf');

Acl = ss.A+ss.B*K;
M = sdpvar(n, n);
delta = 1e-2;
cons = [M(:)- Acl(:) >=0 ; M(:) + Acl(:)>= 0; sum(M, 2) <= lam_test];
opts = sdpsettings;
objective = norm(K, 'inf')*nu/(1-lam_test);

sol = optimize(cons, objective, opts);




%% recovery
Mr = value(M);
Kr = value(K);
Aclr = value(Acl);



Acl_sign = cell(2^m);

Acl_sign{1}= ss.A + ss.B*(eye(2)+diag([1, 1])*delta)*Kr;
Acl_sign{2}= ss.A + ss.B*(eye(2)+diag([1, -1])*delta)*Kr;
Acl_sign{3}= ss.A + ss.B*(eye(2)+diag([-1, 1])*delta)*Kr;
Acl_sign{4}= ss.A + ss.B*(eye(2)+diag([-1, -1])*delta)*Kr;

valid = zeros(2^m,1);
for i=1:2^m
    valid(i) = all(sum(abs(Acl_sign{1}), 2)<= lam_test);
end


if sol.problem==0
    fprintf('feasible with gain %0.4f, want gain<1\n', value(objective))
    if valid
        fprintf('controller is valid\n')
    else
        printf('controller is not valid\n')
    end
else
    fprintf('infeasible\n')
end

% valid_all = all(valid)

% valid = cellfun(@(a) all(sum(abs(a), 2)<lam_test), Acl_sign)


% gainr = value(objective)

%% simulate a trajectory
