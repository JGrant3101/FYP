%% Defining the Hessian for the suspension system of the car w.r.t the x vector NEW.
function J = Suspension_dxdx(x, p)

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
% Mean in log normal function
mew = p(9, :);
% Standard deviation in log normal function
lamda = p(10, :);
% Scaling for Inverse Gaussian distribution
scaling = p(11, :);

% Calculating ride height, this is the static ride height plus the vertical
% displacements of both the sprung and unsprung mass.
h = H + x(1, :) + x(2, :);

% Initialising our Jacobian
J = zeros(4, 4, 4, numel(Zs));

% Populating the non-zero terms with values
J(3, 1, 1, :) = ((-scaling .* exp((-(log(h) - mew).^2) ./ (2 * lamda.^2)) ./ (h.^3 .* Ms .* lamda .* sqrt(2*pi)))) .* ((((mew - log(h)) ./ lamda.^2) - 2)...
    .* (((mew - log(h)) ./ lamda.^2) - 1) - 1./lamda.^2);
J(3, 1, 2, :) = ((-scaling .* exp((-(log(h) - mew).^2) ./ (2 * lamda.^2)) ./ (h.^3 .* Ms .* lamda .* sqrt(2*pi)))) .* ((((mew - log(h)) ./ lamda.^2) - 2)...
    .* (((mew - log(h)) ./ lamda.^2) - 1) - 1./lamda.^2);
J(3, 2, 1, :) = ((-scaling .* exp((-(log(h) - mew).^2) ./ (2 * lamda.^2)) ./ (h.^3 .* Ms .* lamda .* sqrt(2*pi)))) .* ((((mew - log(h)) ./ lamda.^2) - 2)...
    .* (((mew - log(h)) ./ lamda.^2) - 1) - 1./lamda.^2);
J(3, 2, 2, :) = ((-scaling .* exp((-(log(h) - mew).^2) ./ (2 * lamda.^2)) ./ (h.^3 .* Ms .* lamda .* sqrt(2*pi)))) .* ((((mew - log(h)) ./ lamda.^2) - 2)...
    .* (((mew - log(h)) ./ lamda.^2) - 1) - 1./lamda.^2);

end