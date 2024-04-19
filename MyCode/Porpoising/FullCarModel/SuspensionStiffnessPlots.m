ksr = linspace(9*10^4, 2.5*10^5, 1600)';
ksf = ((2.7 * 10^5) * (1.97 * 10^5 - 0.08931 * ksr)) ./ (0.46 * 10^5 + 1.08931 * ksr);

A = ksf ./ ksr;

figure
tiledlayout(1, 2)
nexttile
plot(ksf, ksr)
xlabel('ksf')
ylabel('ksr')
grid minor

nexttile
plot(ksr, A)
xlabel('ksr')
ylabel('A')
grid minor

ksrfinal = interp1(A, ksr, 1.1);
ksffinal = interp1(ksr, ksf, ksrfinal);

