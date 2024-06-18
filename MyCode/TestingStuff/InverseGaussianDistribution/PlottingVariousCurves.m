h = linspace(0, 0.1, 1001)';
original = DWFFloorNew(h, 0.0001, 2.4, 0.31*(500/9)^2); % Original scaling value was 956.7901

lamdas = linspace(1, 3, 41)';
results = zeros(length(h), length(lamdas));
for i = 1:length(lamdas)
    results(:, i) = DWFFloorNew(h, 0.0001, lamdas(i), 0.31*(500/9)^2);
end

figure
hold on
plot(h, original, 'color', 'black')
c1 = [0, 0.447, 0.741];
c2 = [0.85, 0.325, 0.098];
c = [linspace(c1(1),c2(1),length(lamdas))', linspace(c1(2),c2(2),length(lamdas))', linspace(c1(3),c2(3),length(lamdas))'];
for i = 1:length(lamdas)
    plot(h, results(:, i), 'color', c(i, :))
end
legendlabels = {};
legendlabels{1, 1} = "Original";
for i = 1:length(lamdas) 
    legendlabels{i+1, 1} = "lamda = " + num2str(lamdas(i));
end
legend(legendlabels)