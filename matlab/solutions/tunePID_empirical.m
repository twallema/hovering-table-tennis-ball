% clear all variables, close all windows
clc, clear all, close all;

% calibrate system parameters
t_d = 0.5; % s; dead time
alpha = 0.19;
beta = 0.63;

% tuned controllers
C = {[57, 0, 0, 100], [52, 17.8, 0, 100], [67, 38.3, 153, 100]};
labels = {'P', 'PI', 'PID'};

% no noise and no filtering
tau = 0.01; % s; measurement RC filter time constant
noise_power = 0; 
noise_sample_time = 0.05;

% initialise experimental parameters
simtime = 15; % experiment time
z_0 = [0.2, 0]; % initial position and speed
z_sp = 0.4; % desired setpoint (height)
steptime = 5; % time of setpoint change

% select model to use
model = 'ODEModelPID.slx';

% ODE model parameters
rho = 1.225; % kg.m^-3 
A = 4*3.14*0.02^2; % m^2
m = 0.00283; % kg
g = 9.81; % m.s^-2
C_D = 0.47; % -
u_eq = (2*m*g/(rho*A*C_D*alpha^2))^(1/(2*beta)); % eq. fan speed

figure(1)
linestyles = {'-', '-.', ':'};
for i = 1:length(C)
    % unpack controller parameters
    P = C{i}(1);
    I = C{i}(2);
    D = C{i}(3);
    F = C{i}(4);
    % simulate
    load_system(model)
    out = sim(model);
    % extract
    t = out.simout(:,1);
    h = out.simout(:,7);
    % visualise height
    plot(t, h, 'k', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Displayname', labels{i})
    hold on
end
% visualise setpoint
plot(t, out.simout(:,2), 'color', [0.7,0,0,1], 'LineWidth', 2, 'Displayname', 'Setpoint')
hold off
ylim([0,0.65])
xlabel('Time (s)')
ylabel('Height (m)')
legend('Location', 'southeast')
set(gcf,'units','inches','position',[0,0,8.3,0.20*11.7])
% save result
saveas(gcf,'tuned_PID_empirical.pdf')