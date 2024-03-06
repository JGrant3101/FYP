x(1, :) = -0.015;
x(2, :) = -0.005;
x(3, :) = -0.01;
x(4, :) = -0.015;

Ms = 180;
p(2, :) = 50;
Ks = 0.9*10^5;
p(4, :) = 3400;
p(5, :) = 0.9*2.7 * 10^5;
p(6, :) = 0.1;
p(7, :) = ((250.*10^3)./(60*60));
p(8, :) = 0.365;
mew = 0.0001;
lamda = 2.4;
scaling = 0.31*(500/9)^2;

h = p(6, :) + x(1, :) + x(2, :);

origval = -(Ks./Ms) - (scaling./(Ms .* h.^2 .* lamda .* sqrt(2.*pi))) * exp((-(log(h) - mew).^2) ./ (2 * lamda.^2)) * ((mew - log(h))./(lamda.^2) - 1);

hnew = h * 1.00001;

step = h * 0.00001;

newval = -(Ks./Ms) - (scaling./(Ms .* hnew.^2 .* lamda .* sqrt(2.*pi))) * exp((-(log(hnew) - mew).^2) ./ (2 * lamda.^2)) * ((mew - log(hnew))./(lamda.^2) - 1);

deriv = (newval - origval)./step;