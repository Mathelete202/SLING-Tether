%{
Title: Basic Tension
Author: Josiah Kasper
Date: 3/25/24
Description: Models the tension throughout our tether given no payload or
counterweight mechanism.
%}

% declare properties
rho = 13; % kg/m
omega = 2; % rad/s

tension = @(x) rho*omega^2.*x.^2/2;

x = linspace(-100000,100000,1e7);
plot(x/1000,tension(x)/1000)
title("Tension vs. Position")
xlabel("Tether position from C.O.M (km)")
ylabel("Tension (kN)")
grid on
grid minor