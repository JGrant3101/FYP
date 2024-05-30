clear all
close all
clc
k = 3; % Values For ISO Road A-B Roughness Classification, from 3 to 9
N  = 15000; %  Number of data points
L  = 1500;  % Length Of Road Profile (m)
B  = L/N ; % Sampling Interval (m)
dn = 1/L;  % Frequency Band
n0 = 0.1;        % Spatial Frequency (cycles/m)
n  = dn : dn : N*dn; % Spatial Frequency Band
phi =  2*pi*rand(size(n)); % Random Phase Angle
Amp1 = sqrt(dn)*(2^k)*(1e-3)*(n0./n); % Amplitude for Road  Class A-B
x = 0:B:L-B; % Abscissa Variable from 0 to L
hx = zeros(size(x));
for i=1:length(x)
    hx(i) = sum(Amp1.*cos(2*pi*n*x(i)+ phi));
end

Differences = diff(hx);
    
TotalVerticalMove = sum(abs(Differences));
Roadpoints = linspace(0, 1500, 15000);

figure
plot(Roadpoints, hx)

% Filtering only the frequencies of intereset
hxfiltered = bandpass(hx, [0.5, 15], 10);
hxfiltered = hxfiltered(6:end-5);

figure
plot(Roadpoints(6:end-5), hxfiltered)

figure
[q , C] = psd_1D(hxfiltered, B, 'x');  % B is Sampling Interval (m); for the case that I have explained above it will be 250/45000 = 5.55e-3
lambda = (2*pi) ./ q; % wavelengths
f = q / (2*pi); % spatial frequency
PSD = 2 * pi * C; % power spectrum
loglog(f,PSD)