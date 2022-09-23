%% VALIDATION TASK 1

clear all
close all

% Inputs
wind_speed = 2.5;
density = 1.225;
wind_direction = 270;
turbine_centres = [10,0,0.3];
yaw_angles = 0;
diameters = 0.416;
power_curve = csvread('WT_powercurve.csv',1,0);
exp_data = csvread('WT_NeutralABL_exp.csv',1,0);
location = exp_data(301:357,1:3)/1000;
location = exp_data(301:357,1:3)/1000;

% Function floris takes in all inputs and outputs normalised speed
[power, speed] = floris_v2(wind_speed, density, wind_direction, turbine_centres, yaw_angles, diameters, power_curve, location);


% Graph customization
% Compares experimental data and validation data at x5 diameter distance
plot(exp_data(301:357,2)/1000,exp_data(301:357,7),'bluex');
hold on
plot(exp_data(301:357,2)/1000,speed/wind_speed,'black-');
xlabel("Y locations (mm)")
ylabel("Normalized speed predicted by model")
title("Comparison of FLORIS model to wind tunnel data for x5 diameters downstream")
legend("Experimental Data", "Validation Data")

%% VALIDATION TASK 2

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
power_curve = readmatrix('V80_powercurve_1.csv');
location = [0 0 0];

% Function floris takes in all inputs and outputs normalised wind power
[power,speed] = floris_v2(wind_speed, density, wind_direction, turbine_centres, yaw_angles, diameters, power_curve, location);
normalised_power = power/power(1);

% Graph customization
% Plots validation data with FLORIS data for comparison
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