%% Defining the suspension system of the car, this also includes the tyre stiffness.
function [rhs, VarOfInterest] = Suspension(t, x, p)
% Firslty defining the constant terms
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
vCar = p(7, :);
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
h = H + x(1, :) + x(2, :);
% Calculating the downforce from the floor
DWFFloorValue = DWFFloor(h, mew, lamda, scaling);

% Running calculations for vertical displacement first and second derivatives in
% time.
rhs = zeros(size(x));
rhs(1, :) = x(3, :);
rhs(2, :) = x(4, :);
rhs(3, :) = -(Ks./Ms) .* x(1, :) + (Ks./Ms) .* x(2, :) - (Cs./Ms) .* x(3, :) + (Cs./Ms) .* x(4, :) - (DWFUpper + DWFFloorValue)./Ms; % This also has the variable downforce component dependent on Vcar and ride height which will be added
rhs(4, :) = (Ks./Mu) .* x(1, :) - ((Ks./Mu) + (Kt./Mu)) * x(2, :) + (Cs./Mu) * x(3, :) - (Cs./Mu) .* x(4, :);

% Returning values of interest
VarOfInterest(1, 1) = DWFUpper;
VarOfInterest(2, 1) = DWFFloorValue;
VarOfInterest(3, 1) = h;
end


