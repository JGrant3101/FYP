% Optimal damping calculator
ms = 600;
ksf = 1.5 * 2 * 10^5;
ktf = 2 * 2.7 * 10^5;
ksr = 1.5 * 1.81 * 10^5;
ktr = 2 * 0.9 * 2.7 * 10^5;

coptf = sqrt((0.45 * ms * ksf) / 2) * sqrt((ktf + 2*ksf) / ktf);
coptr = sqrt((0.55 * ms * ksr) / 2) * sqrt((ktr + 2*ksr) / ktr);

disp(coptf)
disp(coptr)