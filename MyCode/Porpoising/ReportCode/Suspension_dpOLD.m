%% Defining the Jacobian for the suspension system of the car w.r.t the p vector.
function J = Suspension_dp(x, p)
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
J = zeros(4, 11, numel(Zs));

% Populating the non-zero terms with values
J(3, 1, :) = (1 ./ Ms.^2) .* (Ks .* Zs - Ks .* Zu + Cs .* Phis - Cs .* Phiu + ...
    ((scaling .* sqrt(lamda)) ./ (sqrt(2 * pi) .* h)) .* exp(-(lamda.*(h-mew).^2)./(2.*mew.^2.*h)) + A .* vCar.^2);
J(3, 3, :) = -(Zs ./ Ms) + (Zu ./ Ms);
J(3, 4, :) = -(Phis ./ Ms) + (Phiu ./ Ms);
J(3, 6, :) =  - ((scaling .* sqrt(lamda)) ./ (Ms .* h.^2 .* sqrt(2 .* pi))) .* exp(-(lamda.*(h-mew).^2)./(2.*mew.^2.*h)) .* ...
    (((lamda .* (mew.^2 - h.^2))./(2 .* mew.^2 .* h)) - 1);
J(3, 7, :) = (-2 * A .* vCar) ./ Ms;
J(3, 8, :) = (-vCar.^2) ./ Ms;
J(3, 9, :) = ((-scaling .* lamda.^(3/2) .* (h - mew)) ./ (Ms .* h .* mew.^3 .* sqrt(2 * pi))) .* exp(-(lamda.*(h-mew).^2)./(2.*mew.^2.*h));
J(3, 10, :) = (-scaling ./ (Ms .* h .* sqrt(2 * pi))) .* exp(-(lamda.*(h-mew).^2)./(2.*mew.^2.*h)) .* ...
    (0.5 .* lamda.^(-1/2) - ((sqrt(lamda) .* (h - mew).^2) ./ (2 .* mew.^2 .* h)));
J(3, 11, :) = (-sqrt(lamda) ./ (Ms .* h .* sqrt(2 * pi))) .* exp(-(lamda.*(h-mew).^2)./(2.*mew.^2.*h));

J(4, 2, :) = (1 ./ Mu.^2) .* (-Ks .* Zs + (Ks + Kt) .* Zu - Cs .* Phis + Cs .* Phiu);
J(4, 3, :) = (Zs ./ Mu) - (Zu ./ Mu);
J(4, 4, :) = (Phis ./ Mu) - (Phiu ./ Mu);
J(4, 5, :) = -Zu ./ Mu;
end