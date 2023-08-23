rng(60, 'twister')
% n = 4;
% m = 3;

n = 3;
m = 2;
ss = drss(n, n, m);
% A0 = ss.A + 0.1*eye(n);

SS = 0;
if SS
    ss.A = ss.A*1.4;
    rho = 0.9;
    Kr = [0.528325599264699,-0.819701619416720,-0.517846738580011;-0.851747555414398,0.274528511357424,-0.366597621754505];
    vr = [1, 1, 1];
else
    rho = 0.45;
    ss.A = ss.A*2;
    Kr = [0.283128625069581,-0.623606351013501,-0.284867146189322;0.367523575681024,-1.098900665056812,-2.101179773647827];
    vr = [4.983219006901710;17.691138545433365;8.986643837954036];
end


% Mr = [0.518777602755663,0.652086755295695,0.551198756870782;0.652086755295695,1.016337427112571,0.492662365958030;0.551198756870782,0.492662365958030,0.552589534924221];
% vr = [1.732063114922140;2.171086548366296;1.606450657753033];



%% sample the trajectory
x0 = ones(n, 1);
Tsim = 30;
X = zeros(n, Tsim+1);
U = zeros(m, Tsim);
L = zeros(1, Tsim);

lq = @(u) LogQuant(u, rho);
% lq = @(u) u;

X(:, 1) = x0;
x_curr = x0;
L(1) = norm(x0./vr, 'inf');

for i = 1:Tsim
    %perform the quantization
    u_nominal = Kr*x_curr;
    u_curr = arrayfun(lq, u_nominal);

    x_next = ss.A*x_curr + ss.B*u_curr;
    
    %store the data
    U(:, i) = u_curr;
    X(:, i+1) = x_next;

    L(i+1) = norm(x_next ./ vr, 'inf');
    x_curr = x_next;
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