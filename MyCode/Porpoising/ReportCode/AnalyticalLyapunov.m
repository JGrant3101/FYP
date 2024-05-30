% Trying to get analytic expression for Lyapunov coefficient
clear all
close all
clc
syms Zs Zu Phis Phiu Ms Mu Ks Cs Kt H vCar A mew lamda scaling

% output = Lyapunov(Zs, Zu, Phis, Phiu, Ms, Mu, Ks, Cs, Kt, H, vCar, A, mew, lamda, scaling);

% function output = Lyapunov(Zs, Zu, Phis, Phiu, Ms, Mu, Ks, Cs, Kt, H, vCar, A, mew, lamda, scaling)
x = [Zs; Zu; Phis; Phiu];
p = [Ms; Mu; Ks; Cs; Kt; H; vCar; A; mew; lamda; scaling];

% Find the jacobian
A = Suspension_dx_symbolic(x, p);

% Find the eigenvalues of the Jacobian
[X, D] = eig(A);

% Find the smallet eigenvalue
[m, idx] = min(abs(real(diag(D))));
v = X(:, idx);
om = imag(D(idx, idx));
vb = conj(v);

% Run a check
if m>1e-6
    y = NaN;
    return
end

% Find the eigenvalues of the transpose of the Jacobian
[X, D] = eig(A');

% Find the smallest eigenvalue
[m, idx] = min(abs(real(diag(D))));

w = X(:, idx);

if om*imag(D(idx, idx))>0
    w = conj(w);
end
w = w/conj(w'*v);

B = Suspension_dxdx(x, p);
B1 = zeros(numel(x), 1);
B3 = zeros(numel(x), 1);

for i = 1:numel(x)
    Bmat = reshape(B(i, :, :), [numel(x), numel(x)]);
    B1(i) = v.'*Bmat*v;
    B3(i) = v.'*Bmat*vb;
end


t1 = (2*sqrt(-1)*om*eye(numel(x))-A)\B1;
t2 = A\B3;
B2 = zeros(numel(x), 1);
B4 = zeros(numel(x), 1);

for i = 1:numel(x)
    Bmat = reshape(B(i, :, :), [numel(x), numel(x)]);
    B2(i) = vb.'*Bmat*t1;
    B4(i) = v.'*Bmat*t2;
end

y = real(w'*B2-2*w'*B4)/2/om;
% end