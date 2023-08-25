
rho = 0.4;

delta = (1-rho)/(1+rho);

U = 8;
N = 10000;
u = linspace(-U, U, N);
y = zeros(size(u));
for i = 1:length(u)
    y(i) = LogQuant(u(i), rho);
end

ub= (1+delta)*u;
lb= (1-delta)*u;

figure(1)
tiledlayout(2, 1)
ax1 = nexttile;


c = linspecer(3);


fsa = 14;
hold on
plot(u, ub, 'k', 'LineWidth', 2);
plot(u, lb, 'k', 'LineWidth', 2);
plot(u, u, 'color', c(1, :), 'LineWidth', 2);
plot(u, y, 'color', c(2, :), 'LineWidth', 2);

xlabel('$u$', 'Interpreter','latex', 'fontsize', fsa);
ylabel('$g_{\rho}(u)$', 'Interpreter', 'Latex', 'fontsize', fsa)
title('Logarithmic Quantizer', 'fontsize', 16)
legend({'sector bounds', '', 'nominal', 'quantized'}, 'location', 'southeast', 'fontsize', 10)



ax2 = nexttile;
hold on
plot(u, (y-u)./u, 'color', c(3, :), 'LineWidth', 2)
plot(xlim, -delta*[1, 1], 'k', 'LineWidth', 3);
plot(xlim, delta*[1, 1], 'k', 'LineWidth', 3);
legend('error', '\pm \delta', '', 'location', 'northeast', 'fontsize', 10)
xlabel('$u$', 'Interpreter','latex', 'fontsize', fsa);
ylabel('$(g_{\rho}(u)-u)/u$', 'Interpreter', 'Latex', 'fontsize', fsa)
title('Bounded Multiplicative Error', 'fontsize', 16)
