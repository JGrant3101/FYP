syms zs zu phis phiu ms mu ks cs kt H vCar A mew std scaling
F = -(ks/ms) * zs + (ks/ms) * zu - (cs/ms) * phis + (cs/ms) * phiu - ...
    (1/ms) * (A * vCar^2 + (scaling/((H+zs+zu)*std*sqrt(2*pi)) ...
    * exp(-(log(H+zs+zu) - mew))));