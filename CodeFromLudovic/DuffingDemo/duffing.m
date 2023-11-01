%% Definition of right-hand side of the Duffing equation
function rhs = duffing(t, x, p)
% The right-hand side of the ordinary differential equation
%
% $$ \ddot x + 2\xi\dot x + \omega^2x + \beta x^3 = \Gamma\sin(\Omega t) $$
%
% where $p = [\xi, \omega, \beta, \Gamma, \Omega]$

% By default COCO expects the vector field to be vectorised

xi = p(1, :);
omega = p(2, :);
beta = p(3, :);
Gamma = p(4, :);
Omega = p(5, :);

% The right-hand-side of the differential equations is scaled such that
% $\dot x = Tf(x)$ where $T$ is the period of oscillation. Hence the period
% for COCO is simply one.
T = 2*pi./Omega;

rhs = zeros(size(x));
rhs(1, :) = T.*x(2, :);
rhs(2, :) = T.*(Gamma.*sin(2*pi*t) - 2*xi.*x(2, :) - omega.^2.*x(1, :) - beta.*x(1, :).^3);

end