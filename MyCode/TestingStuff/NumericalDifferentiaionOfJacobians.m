% Numerically differentiating to check for errors
%% Start by initialising input values
clear all
x(1, :) = -0.015;
x(2, :) = -0.005;
x(3, :) = -0.01;
x(4, :) = -0.015;

p(1, :) = 180;
p(2, :) = 50;
p(3, :) = 0.9*10^5;
p(4, :) = 3400;
p(5, :) = 0.9*2.7 * 10^5;
p(6, :) = 0.1;
p(7, :) = ((250.*10^3)./(60*60));
p(8, :) = 0.365;
p(9, :) = 0.0001;
p(10, :) = 2.4;
p(11, :) = 0.31*(500/9)^2;

% Calculating downforce from upper aero elements
% DWFUpper = A.* vCar.^2;
% Calculating ride height, this is the static ride height plus the vertical
% displacements of both the sprung and unsprung mass.
% h = H + Zs + Zu;
% Calculating the downforce from the floor
% DWFFloorValue = DWFFloor(h, mew, lamda, scaling);
% 
% rhs(1, :) = Phis;
% rhs(2, :) = Phiu;
% rhs(3, :) = -(Ks./Ms) .* Zs + (Ks./Ms) .* Zu - (Cs./Ms) .* Phis + (Cs./Ms) .* Phiu - (DWFUpper + DWFFloorValue)./Ms;
% rhs(4, :) = (Ks./Mu) .* Zs - ((Ks./Mu) + (Kt./Mu)) .* Zu + (Cs./Mu) .* Phis - (Cs./Mu) .* Phiu;

%% Jx
J = zeros(4, 4, 1);
for i = 1:4
    for j = 1:4
        if i==1
            % Find initial value of function
            x0 = x(3, :);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* x(j, :);
            temp = x;
            temp(j, :) = 1.000000001 .* x(j, :);
            % Find the new estimate of the function
            x1 = temp(3, :);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        elseif i==2
            % Find initial value of function
            x0 = x(4, :);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* x(j, :);
            temp = x;
            temp(j, :) = 1.000000001 .* x(j, :);
            % Find the new estimate of the function
            x1 = temp(4, :);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        elseif i==3
            DWFUpper = p(8,:).* p(7,:).^2;
            h = p(6,:) + x(1,:) + x(2,:);
            DWFFloorValue = DWFFloorNew(h, p(9,:), p(10,:), p(11,:));
            % Find initial value of function
            x0 = -(p(3,:)./p(1,:)) .* x(1,:) + (p(3,:)./p(1,:)) .* x(2,:) - (p(4,:)./p(1,:)) .* x(3,:) + (p(4,:)./p(1,:)) .* x(4,:) - (DWFUpper + DWFFloorValue)./p(1,:);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* x(j, :);
            temp = x;
            temp(j, :) = 1.000000001 .* x(j, :);
            htemp = p(6,:) + temp(1,:) + temp(2,:);
            DWFFloorValue = DWFFloorNew(htemp, p(9,:), p(10,:), p(11,:));
            % Find the new estimate of the function
            x1 = -(p(3,:)./p(1,:)) .* temp(1,:) + (p(3,:)./p(1,:)) .* temp(2,:) - (p(4,:)./p(1,:)) .* temp(3,:) + (p(4,:)./p(1,:)) .* temp(4,:) - (DWFUpper + DWFFloorValue)./p(1,:);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        elseif i==4
            % Find initial value of function
            x0 = (p(3,:)./p(2,:)) .* x(1,:) - ((p(3,:)./p(2,:)) + (p(5,:)./p(2,:))) .* x(2,:) + (p(4,:)./p(2,:)) .* x(3,:) - (p(4,:)./p(2,:)) .* x(4,:);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* x(j, :);
            temp = x;
            temp(j, :) = 1.000000001 .* x(j, :);
            % Find the new estimate of the function
            x1 = (p(3,:)./p(2,:)) .* temp(1,:) - ((p(3,:)./p(2,:)) + (p(5,:)./p(2,:))) .* temp(2,:) + (p(4,:)./p(2,:)) .* temp(3,:) - (p(4,:)./p(2,:)) .* temp(4,:);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        end
    end
end

%% Jp
J = zeros(4, 11, 1);
for i = 1:4
    for j = 1:11
        if i==1
            % Find initial value of function
            x0 = x(3, :);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* p(j, :);
            temp = p;
            temp(j, :) = 1.000000001 .* p(j, :);
            % Find the new estimate of the function
            x1 = x(3, :);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        elseif i==2
            % Find initial value of function
            x0 = x(4, :);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* p(j, :);
            temp = p;
            temp(j, :) = 1.000000001 .* p(j, :);
            % Find the new estimate of the function
            x1 = x(4, :);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        elseif i==3
            DWFUpper = p(8,:).* p(7,:).^2;
            htemp = p(6,:) + x(1,:) + x(2,:);
            DWFFloorValue = DWFFloor(htemp, p(9,:), p(10,:), p(11,:));
            % Find initial value of function
            x0 = -(p(3,:)./p(1,:)) .* x(1,:) + (p(3,:)./p(1,:)) .* x(2,:) - (p(4,:)./p(1,:)) .* x(3,:) + (p(4,:)./p(1,:)) .* x(4,:) - (DWFUpper + DWFFloorValue)./p(1,:);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* p(j, :);
            temp = p;
            temp(j, :) = 1.000000001 .* p(j, :);
            % Find the new estimate of the function
            DWFUpper = temp(8,:).* temp(7,:).^2;
            htemp = temp(6,:) + x(1,:) + x(2,:);
            DWFFloorValue = DWFFloor(htemp, temp(9,:), temp(10,:), temp(11,:));
            x1 = -(temp(3,:)./temp(1,:)) .* x(1,:) + (temp(3,:)./temp(1,:)) .* x(2,:) - (temp(4,:)./temp(1,:)) .* x(3,:) + (temp(4,:)./temp(1,:)) .* x(4,:) - (DWFUpper + DWFFloorValue)./temp(1,:);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        elseif i==4
            % Find initial value of function
            x0 = (p(3,:)./p(2,:)) .* x(1,:) - ((p(3,:)./p(2,:)) + (p(5,:)./p(2,:))) .* x(2,:) + (p(4,:)./p(2,:)) .* x(3,:) - (p(4,:)./p(2,:)) .* x(4,:);
            % Make the change to the relevant term for this loop
            h = 0.000000001 .* p(j, :);
            temp = p;
            temp(j, :) = 1.000000001 .* p(j, :);
            % Find the new estimate of the function
            x1 = (temp(3,:)./temp(2,:)) .* x(1,:) - ((temp(3,:)./temp(2,:)) + (temp(5,:)./temp(2,:))) .* x(2,:) + (temp(4,:)./temp(2,:)) .* x(3,:) - (temp(4,:)./temp(2,:)) .* x(4,:);
            % Calculate the derivative
            deriv = (x1 - x0)./h;
            % Put the value in the Jacobian
            J(i, j) = deriv;
        end
    end
end