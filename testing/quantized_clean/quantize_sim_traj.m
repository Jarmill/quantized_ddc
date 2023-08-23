rng(30, 'twister')
% n = 4;
% m = 3;

n = 3;
m = 2;
ss = drss(n, n, m);
% A0 = ss.A + 0.1*eye(n);
ss.A = ss.A*1.5;

q = 0.1665;
rho = (1+q)/(1-q);

Mr = [0.518777602755663,0.652086755295695,0.551198756870782;0.652086755295695,1.016337427112571,0.492662365958030;0.551198756870782,0.492662365958030,0.552589534924221];
vr = [1.732063114922140;2.171086548366296;1.606450657753033];
Kr = [0.274453503818638,0.330651120472518,-0.513421004212335;0.059223975414846,0.267926983857520,-0.077982352955193];


%% sample the trajectory
x0 = ones(n, 1);
Tsim = 30;
X = zeros(n, Tsim+1);
U = zeros(m, Tsim);
L = zeros(1, Tsim);

lq = @(u) LogQuant(u, q);

X(:, 1) = x0;
x_curr = x0;
L(1) = norm(x0./vr, 'inf');

for i = 1:Tsim
    %perform the quantization
    u_nominal = Kr*x_curr;
    u_curr = arrayfun(lq, u_nominal);

    x_next = ss.A * x_curr + ss.B*u_curr;
    
    %store the data
    U(:, i) = u_curr;
    X(:, i+1) = x_next;

    L(i+1) = norm(x_next ./ vr, 'inf');
end

%% plot the results

figure(1)
plot(X')
xlabel('t')
ylabel('x')

figure(2)
plot(L')
xlabel('t')
ylabel('norminf(x./v)')