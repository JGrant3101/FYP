% Fitting quadratic regression to family of Hopf points
vCar = coco_bd_col(HBbd{1}, 'vCar')';
otherparam = coco_bd_col(HBbd{1}, param4HB)';
otherparamnormalised = otherparam / Inputs(Index4HB, 1);
p = polyfit(otherparamnormalised, vCar, 4);

vCarPredictions = p(1) * otherparamnormalised.^4 + p(2) * otherparamnormalised.^3 + p(3) * otherparamnormalised.^2 + p(4) * otherparamnormalised + p(5);

figure
scatter(otherparamnormalised, vCar)
hold on
plot(otherparamnormalised, vCarPredictions, 'LineWidth', 3)
xlabel('Normalised Ks')
ylabel('vCar (kph)')
legend('COCO', 'Polynomial')
title('vCar vs normalised Ks')

percentdiff = ((vCarPredictions - vCar) ./ vCar) * 100;
avgpercentdiff = mean(percentdiff);
disp(avgpercentdiff)

Derivs = 4 * p(1) * otherparamnormalised.^3 + 3 * p(2) * otherparamnormalised.^2 + 2 * p(3) * otherparamnormalised + p(4);

vCarfilt = vCar(vCar > 200);
Derivs = Derivs(vCar > 200);

avgDeriv = mean(Derivs);
disp(avgDeriv)