%% Defining the suspension system of the car, this also includes the tyre stiffness.
function rhs = SuspensionWithTime(t, x, p)

% Firstly define the variable terms
% Front unsprung mass vertical displacement
ZuF = x(1, :);
% Rear unsprung mass vertical displacement
ZuR = x(2, :);
% Sprung mass vertical displacement
Zs = x(3, :);
% Sprung mass angle of rotation
THETAs = x(4, :);

% Front unsprung mass vertical velocity
PhiuF = x(5, :);
% Front unsprung mass vertical velocity
PhiuR = x(6, :);
% Front unsprung mass vertical velocity
Phis = x(7, :);
% Front unsprung mass vertical velocity
PhiTHETAs = x(8, :);

% Then defining the constant terms
% Unsprung front mass
MuF = p(1, :);
% Unsprung rear mass
MuR = p(2, :);
% Sprung mass
Ms = p(3, :);
% Sprung mass moment of inertia
Is = p(4, :);
% Distance from front axle to CoG
LF = p(5, :);
% Distance from rear axle to CoG
LR = p(6, :);
% Front tyre vertical stiffness
KtF = p(7, :);
% Front suspension stiffness
KsF = p(8, :);
% Rear tyre vertical stiffness
KtR = p(9, :);
% Rear suspension stiffness
KsR = p(10, :);
% Front damping 
CsF = p(11, :);
% Rear damping
CsR = p(12, :);
% Static front ride height
StaticFRH = p(13, :);
% Static rear ride height
StaticRRH = p(14, :);
% Car speed
vCar = ((p(15, :).*10^3)./(60*60));
% Front wing aero multiplier
AF = p(16, :);
% Rear wing aero multiplier
AR = p(17, :);
% Front downforce from floor gradient
m = p(18, :);
% Front downforce from floor constant
DWFFloorFConstant = p(41, :);

% Downforce from floor polynomial terms
p00 = p(19, :);
p10 = p(20, :);
p01 = p(21, :);
p20 = p(22, :);
p11 = p(23, :);
p02 = p(24, :);
p30 = p(25, :);
p21 = p(26, :);
p12 = p(27, :);
p03 = p(28, :);
%p40 = p(29, :);
p31 = p(29, :);
p22 = p(30, :);
p13 = p(31, :);
p04 = p(32, :);
%p50 = p(34, :);
%p41 = p(35, :);
p32 = p(33, :);
p23 = p(34, :);
p14 = p(35, :);
p05 = p(36, :);
avgRHmean = p(37, :);
avgRHstd = p(38, :);
vCarmean = p(39, :);
vCarstd = p(40, :);

% Floor road interaction parameters
e = p(42, 1);

% Calculating downforce values from wings
DWFFrontWing = AF .* vCar.^2;
DWFRearWing = AR .* vCar.^2;

% Calculating dynamic ride heights
hF = StaticFRH + Zs + ZuF - LF .* THETAs;

% Defining the floor geometry
FloorLength = 3;
FloorIndices = 30001;

% Create an array of floor distances based on the floor length
FloorArray = linspace(0, FloorLength, FloorIndices)';

% Defining the x and y for the low point of the floor
Xlow = 1.7;
Ylow = StaticRRH - 0.03;

% Calculating the static rake angle
StaticRake = atan((StaticRRH - StaticFRH) / (LF + LR));

% Will form a cubic defining the floor profile
temp1 = zeros([4, 4]);
temp2 = zeros([4, 1]);

temp1(1, 4) = 1;
temp1(2, 1) = FloorLength^3;
temp1(2, 2) = FloorLength^2;
temp1(2, 3) = FloorLength;
temp1(2, 4) = 1;
temp1(3, 1) = Xlow^3;
temp1(3, 2) = Xlow^2;
temp1(3, 3) = Xlow;
temp1(3, 4) = 1;
temp1(4, 1) = 3 * Xlow^2;
temp1(4, 2) = 2 * Xlow;
temp1(4, 3) = 1;

temp2(1) = StaticFRH + 0.6 * tan(StaticRake);
temp2(2) = StaticRRH;
temp2(3) = Ylow;
temp2(4) = 0;

% Solving
values = linsolve(temp1, temp2);
a = values(1);
b = values(2);
c = values(3);
d = values(4);

% Defining the floor profile
StaticFloorProfile = a * FloorArray.^3 + b * FloorArray.^2 + c * FloorArray + d;

% Converting this static floor profile into the dynamic floor profile
% First by applying the sprung mass displacment
FloorProfile = StaticFloorProfile + Zs;

% Then by applying the angle of the sprung mass
THETAChange = THETAs - atan((StaticRRH - StaticFRH) / (LF + LR));
% Need to find the distance of each point in the floor from the point of
% rotation to find it's veritcal displacment
FloorDist2Rotation = FloorArray - (LF - 0.6);
% Now finding the veritcal displacement due to the dynamic rotation of the
% floor
FloorVertChangeDue2Rotation = FloorDist2Rotation * tan(THETAChange);
% Applying these vertical changes
FloorProfile = FloorProfile + FloorVertChangeDue2Rotation;
% Now that the dynamic floor profile has been found can compute the average
% distance from the road to the floor
avgRH = mean(FloorProfile);
FloorProfileMin = min(FloorProfile);

% Normalising the vCar and average distance to road values 
avgRHnormalised = (avgRH - avgRHmean) ./ avgRHstd;
vCarnormalised = (p(15, :) - vCarmean) ./ vCarstd;

% Setting variables as x and y 
z = avgRHnormalised;
y = vCarnormalised;

% Predicting the total downforce produced by the floor
TotalDWFFloor = p00 + p10.*z + p01.*y + p20.*z.^2 + p11.*z.*y + p02.*y.^2 + p30.*z.^3 + p21.*z.^2.*y + p12.*z.*y.^2 + p03.*y.^3 ... % + p40.*z.^4 ...
    + p31.*z.^3.*y + p22.*z.^2.*y.^2 + p13.*z.*y.^3 + p04.*y.^4 ... % + p50.*z.^5 + p41.*z.^4.*y ...
    + p32.*z.^3.*y.^2 + p23.*z.^2.*y.^3 + p14.*z.*y.^4 + p05.*y.^5;

% TotalDWFFloor = heaviside(FloorProfileMin) .* TotalDWFFloor;

% Calculating the DWF on the front axle from the floor using the linear
% approximation
% DWFFloorF = DWFFloorFConstant + m .* vCar.^2 .* (StaticFloorProfile(1) - FloorProfile(1));
DWFFloorF = 0.45 * TotalDWFFloor;
% DWFFloorF = (DWFFloorFConstant + m .* ((StaticFloorProfile(1) - FloorProfile(1)) ./ StaticFloorProfile(1))) .* TotalDWFFloor;

% From this obtaining the DWF on the rear axle from the floor
DWFFloorR = TotalDWFFloor - DWFFloorF;

% Therefore calculating the total DWF acting on the front and rear axles
DWFF = DWFFrontWing + DWFFloorF;
DWFR = DWFRearWing + DWFFloorR;
% DWFF = DWFFrontWing;
% DWFR = DWFRearWing;

% Can now run the parameters through the equations of motion
rhs = zeros(size(x));
rhs(1, :) = PhiuF;
rhs(2, :) = PhiuR;
% rhs(3, :) = heaviside(FloorProfileMin) .* Phis + heaviside(-FloorProfileMin) .* -e .* Phis;
rhs(3, :) = Phis;
rhs(4, :) = PhiTHETAs;

rhs(5, :) = (KsF ./ MuF) .* (Zs - ZuF - LF .* THETAs) + (CsF ./ MuF) .* (Phis - PhiuF - LF .* PhiTHETAs) - (KtF ./ MuF) .* ZuF;
rhs(6, :) = (KsR ./ MuR) .* (Zs + LR .* THETAs - ZuR) + (CsR ./ MuR) .* (Phis + LR .* PhiTHETAs - PhiuR) - (KtR ./ MuR) .* ZuR;
% rhs(7, :) = -(KsF ./ Ms) .* (Zs - ZuF - LF .* THETAs) - (CsF ./ Ms) .* (Phis - PhiuF - LF .* PhiTHETAs) - (KsR ./ Ms) .* (Zs + LR .* THETAs - ZuR) ...
%     - (CsR ./ Ms) .* (Phis + LR .* PhiTHETAs - PhiuR) - DWFF ./ Ms - DWFR ./ Ms + heaviside(-FloorProfileMin) .* Phis .* (e + 1);
rhs(7, :) = -(KsF ./ Ms) .* (Zs - ZuF - LF .* THETAs) - (CsF ./ Ms) .* (Phis - PhiuF - LF .* PhiTHETAs) - (KsR ./ Ms) .* (Zs + LR .* THETAs - ZuR) ...
    - (CsR ./ Ms) .* (Phis + LR .* PhiTHETAs - PhiuR) - DWFF ./ Ms - DWFR ./ Ms;
rhs(8, :) = (LF ./ Is) .* (KsF .* (Zs - ZuF - LF .* THETAs) + CsF .* (Phis - PhiuF - LF .* PhiTHETAs) + DWFF) ...
    - (LR ./ Is) .* (KsR .* (Zs + LR .* THETAs - ZuR) + CsR .* (Phis + LR .* PhiTHETAs - PhiuR) + DWFR);

end