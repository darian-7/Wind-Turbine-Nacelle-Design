clear all
close all
clc

horizontal_position = [0; 560; 1120; 1680; 2240; 2800; 3360; 3920; 4480; 5040];
normalised_power_validation = [1; 0.68; 0.66; 0.63; 0.62; 0.61; 0.59; 0.57; 0.56; 0.54];


% Inputs

wind_speed = 8;
density = 1.225;
wind_direction = 180;
diameters = 80*ones(10,1);
turbine_center_x = zeros(10,1);
turbine_center_y = horizontal_position;
turbine_center_z = diameters;
turbine_centres = [turbine_center_x turbine_center_y turbine_center_z];
yaw_angles = zeros(10,1);
power_curve = csvread('V80_powercurve.csv',1,0);
location = [0 0 0];


% Function floris takes in all inputs and outputs normalised wind power

[power, speed] = floris(wind_speed, density, wind_direction, turbine_centres, yaw_angles, diameters, power_curve, location);
normalised_power = power/power(1);


% Graph customization

figure
hold on
box on
grid on
plot(horizontal_position, normalised_power_validation, '-b')
plot(horizontal_position, normalised_power, '--k')
xlabel("Position of turbine centres (m)")
ylabel("Normalised Power (W)")
title("Normalised Power Output of Wind Turbines")
legend("Validation Data", "FLORIS Calculated Data")