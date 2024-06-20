t = tiledlayout(2,2); % Function that allows all four graphs to be plotted simueltaneously

%% WIND FARM LAYOUT

clear all
close all
clc

% Iteration that assigns each wind turbine with three distinct coordinates

% Initialization of variables
n = 1; 
c = 550;
d = 0;

pos = zeros(3,144);

for i = 1:143;
    pos(1,n) = c; % x value
    pos(3,n) = 100; % z value
    n=n+1;
    c=c+550;
    pos(2,n) = d; % y value
    if c>6050;
        d = d+770;
        c = 0;
    end
end

% Correctional value
pos(3,1)=100;
pos(3,144)=100;

% Graph plotting
disp(pos);
A = pos(1,:);
B = pos(2,:);

% Graph customization
% Plots a grid of an aerial view ofthe wind turbine farm
hold on
box on

nexttile
scatter(A,B,25, 'x');
xlabel("x coordinates of length abreast")
ylabel("y coordinates of length downwind")
title("Wind farm layout")
box on;
xlim([-400 6400]);
ylim([-600 8900]);

%% WIND SPEEDS THROUGH THE FARM FOR ONE GIVEN WIND SPEED

clear all
clc

% Initialization of variables
n = 1;
c = 550;
d = 0;

turbine_centres = zeros(144,3);

for i = 1:143;
    turbine_centres(n,1) = c; % x value
    turbine_centres(n,3) = 100; % y value
    n=n+1;
    c=c+550;
    turbine_centres(n,2) = d; % z value 
    if c>6050;
        d = d+770;
        c = 0;
    end
end

% Correctional value
turbine_centres(1,3)=100;
turbine_centres(144,3)=100;

% Inputs for the power curve
wind_speed = 21; % calculated by taking the average median wind speed over 12 months
density = 1.225;
wind_direction = 180;
diameters = 110*ones(144,1);
power_curve = readmatrix('V80_powercurve_2.csv');
yaw_angles = zeros(144,1);

% Initialization of variables
a = 1;
x = 0; 
y = 0;

% Locations behind first turbine
for i = 1:40;
    location(a,1) = x; % x value
    location(a,3) = 100; % z value
    location(a,2) = y; % y value
    a=a+1;
    if x == 0;
        y = y+100;
        x = 0;
    end
end

% Calling the floris function
[power,speed] = floris_v2(wind_speed, density, wind_direction, turbine_centres, yaw_angles, diameters, power_curve, location);

% Graph customization
% Plots how wind speed varies in turbines arranged in a row downwind
nexttile
hold on
box on
grid on
plot(location(:,2), speed, 'blue-x')
title('Wind speeds through farm for one given wind speed')
xlabel('Downwind distance from first turbine (m)')
ylabel('Wind speed (m/s)')

%% EFFECT OF CHANGING ONCOMING FLOW TO POWER OUTPUT

clear all
clc

% Initialization of variables
n = 1;
c = 550;
d = 0;

turbine_centres = zeros(140,3);

for i = 1:139;
    turbine_centres(n,1) = c; % x value
    turbine_centres(n,3) = 100; % z value
    n=n+1;
    c=c+550;
    turbine_centres(n,2) = d; % y value
    if c>6050;
        d = d+770;
        c = 0;
    end
end

% Correctional value
turbine_centres(1,3)=100;
turbine_centres(140,3)=100;

% Inputs
wind_speed = 15;
density = 1.225;

diameters = 110*ones(140,1);

power_curve = csvread('v80_powercurve_2.csv',1,0);

location = [0,0,0];

yaw_angles = zeros(140,1);

% Initialization for wind direction variable
wd = 0; 

for i = 1:360
    wind_direction = wd; % Sets wind direction to variable 'wd'
    [power] = floris(wind_speed, density, wind_direction, turbine_centres, yaw_angles, diameters, power_curve, location);
    total_power = (sum(power));
    p_power(i,:)=[total_power];
    wd=wd+1; % Increments wind direction by 1 iteration to complete 360 degrees
end

% Graph cuztomization
% Plots how the power changes with wind direction (measured in degrees)
nexttile
hold on
box on
grid on
plot(0:359,p_power, '-k')
xlim([0 360])
title("Total Farm Power Output against Wind Direction")
xlabel("Degrees")
ylabel("Power (W)")

%% POWER OUTPUT AT DIFFERENT MONTHS IN THE YEAR

clear all
clc

% Inputs wind speed
wind_speeds = readmatrix('average_wind_speed.csv');

% Initialization
n = 1;
c = 550;
d = 0;

for i = 1:143;
    turbine_centres(n,1) = c; % x value
    turbine_centres(n,3) = 100; % z value
    n=n+1;
    c=c+550;
    turbine_centres(n,2) = d; % y value
    if c>6050;
        d = d+770;
        c = 0;
    end
end

% Correctional value
turbine_centres(1,3)=100;
turbine_centres(144,3)=100;

for i = 1:12 % Iteration for 12 months
    wind_speed = wind_speeds(i); % Selecting wind speed from the matrix
    wind_direction = 180;
    if i == 6 || i == 7
        wind_direction = 270; % Since different months have different wind direction
    end 
    density = 1.225*1.05;
    if i > 3 && i < 7
        density = 1.225;
    elseif i >=7 && i <= 10
        density = 1.225*0.99;
    elseif i > 10
        density = 1.225; % Since different months have different density
    end
    power_curve = csvread('v80_powercurve_2.csv',1,0);
    location = [0 0 0];
    yaw_angles = zeros(144,1);
    diameters = 110*ones(144,1);

% Uses FLORIS function to model the power and wind speed of each turbine at their specified location    
  [power] = floris_v2(wind_speed, density, wind_direction, turbine_centres, yaw_angles, diameters, power_curve, location); 
 
  total_power(1,i) = sum(power,'all')

 
end

% Categorizing different power outputs as per their respective month
m = categorical({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
months = reordercats(m,{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});

% Graph customization
% Plots the predicted monthly power output throughout the year

nexttile
hold on
grid on
box on

bar(months,total_power/1E6)
title('Variation of power output during the year')
xlabel('Month')
ylabel('Monthly power output (MW)')
ylim([0 1200]);
yearly_output = sum(total_power)/(12*1E6);
disp('Average power output per annum (MW) = ')
disp(yearly_output)
yline(yearly_output) 
text(4, 650, "Average power output per annum = 611 MW")