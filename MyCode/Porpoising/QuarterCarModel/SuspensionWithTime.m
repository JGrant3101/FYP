function rhs  = SuspensionWithTime(t, x, p)

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

% Calculating downforce from upper aero elements
DWFUpper = A.* vCar.^2;
% Calculating ride height, this is the static ride height plus the vertical
% displacements of both the sprung and unsprung mass.
h = H + Zs + Zu;
% Calculating the downforce from the floor
DWFFloorValue = DWFFloor(h, mew, lamda, scaling);

% Running calculations for vertical displacement first and second derivatives in
% time.
rhs = zeros(size(x));
rhs(1, :) = Phis;
rhs(2, :) = Phiu;
rhs(3, :) = -(Ks./Ms) .* Zs + (Ks./Ms) .* Zu - (Cs./Ms) .* Phis + (Cs./Ms) .* Phiu - (DWFUpper + DWFFloorValue)./Ms;
rhs(4, :) = (Ks./Mu) .* Zs - ((Ks./Mu) + (Kt./Mu)) .* Zu + (Cs./Mu) .* Phis - (Cs./Mu) .* Phiu;

end