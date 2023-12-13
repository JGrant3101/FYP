%% Defining the Jacobian for the suspension system of the car w.r.t the x vector.
function J = Suspension_dx(x, p)

% Firstly define the variable terms
% Sprung mass vertical displacement
Zs = x(1, :);
% Unsprung mass vertical displacement
Zu = x(2, :);
% Sprung mass vertical velocity
Phis = x(3, :);
% Unsprung mass vertical velocity
Phiu = x(4, :);

% Then defining the constant terms
% Sprung mass
Ms = p(1, :);
% Unsprung mass
Mu = p(2, :);
% Suspension stiffness
Ks = p(3, :);
% Suspension damping
Cs = p(4, :);
% Tyre stiffness
Kt = p(5, :);
% Static ride height
H = p(6, :);
% Car speed
vCar = ((p(7, :).*10^3)./(60*60));
% Upper downforce elements multiplier
A = p(8, :);
% Mew for Inverse Gaussian distribution
mew = p(9, :);
% Lamda for Inverse Gaussian distribution
lamda = p(10, :);
% Scaling for Inverse Gaussian distribution
scaling = p(11, :);

% Calculating ride height, this is the static ride height plus the vertical
% displacements of both the sprung and unsprung mass.
h = H + x(1, :) + x(2, :);

% Initialising our Jacobian
J = zeros(4, 4, numel(Zs));

% Populating the non-zero terms with values
J(1, 3, :) = 1;

J(2, 4, :) = 1;

J(3, 1, :) = -(Ks./Ms) - ((scaling .* sqrt(lamda)) ./ (Ms .* sqrt(2 .* pi))) .* ((lamda .* (h - mew) .* ...
    ((h - mew) - 2 .* h) - 2 .* mew.^2) ./ 2 .* mew.^2 .* h.^2) .* exp(-(lamda.*(h-mew).^2)./(2.*mew.^2.*h));
J(3, 2, :) = (Ks./Ms) - ((scaling .* sqrt(lamda)) ./ (Ms .* sqrt(2 .* pi))) .* ...
    ((lamda .* (h - mew) .* ((h - mew) - 2 .* h) - 2 .* mew.^2) ./ 2 .* mew.^2 .* h.^2) .* exp(-(lamda.*(h-mew).^2)./(2.*mew.^2.*h));
J(3, 3, :) = -Cs ./ Ms;
J(3, 4, :) = Cs ./ Ms;

J(4, 1, :) = Ks ./ Mu;
J(4, 2, :) = -((Ks ./ Mu) + (Kt ./ Mu));
J(4, 3, :) = Cs ./ Mu;
J(4, 4, :) = -Cs ./ Mu;

end