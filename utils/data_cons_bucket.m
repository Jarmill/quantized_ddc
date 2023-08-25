function [C, d] = data_cons_bucket(sim, reduce_face)
%generate inequality constraints 
%norm(Xdelta - AX - Bu)_infty <= epsilon
%
%polytope description:
% C[vec(A); vec(B)] <= d
%extract the data

% sim = struct('X', Xn, 'U', U, 'buckets', buckets);
% sim.Sb =  Sb;

if nargin < 2
    reduce_face=1;
end

m = size(sim.U, 1);
[n, T] = size(sim.X);

Xn = sim.X;
U = sim.U;

%get bounds from the buckets
lb = zeros(n, T);
ub = zeros(n, T);
Nbucket = length(sim.Sb);
for i = 1:Nbucket
    lb(sim.Sb{i}) = sim.buckets(i, 1);
    ub(sim.Sb{i}) = sim.buckets(i, 2);
end

mask_lb_inf = (lb~=-Inf);
mask_ub_inf = (ub~=Inf);
%form expressions
CA0 = kron(Xn', speye(n));
CB0 = kron(U', speye(n));

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
if reduce_face
    [C_red, d_red] = nontrivial_constraints(C, d);
    C = C_red; 
    d = d_red;
end


end