% clear all variables, close all windows
clc, clear all, close all;

% controller parameters
K_P = 300; % proportional
K_I = 10; % integral
K_D = 30; % derivative
K_F = 5; % derivative time constant

% calibrated system parameters
t_d = 0.5; % s; dead time
alpha = 0.1885;
beta = 0.6272;

% initialise experimental parameters
simtime = 40; % experiment time
z_0 = [0, 0]; % initial position and speed
z_sp = 0.4; % desired setpoint (height)
steptime = 5; % time of setpoint change
disturb_time = 20; % time of fan effectiveness manipulation
disturb_magnitude = -7.5; % percentage change in parameter alpha per unit time

% noise and filter parameters
tau = 0.05; % s; measurement RC filter time constant
noise_power = 0; % 0.0008; 
noise_sample_time = 0.05;

% select model to use
model = '../to_students/ODEModelPID_smithpredictor.slx'; % alternative: ../to_students/ODEModelPID_smithpredictor.slx

% ODE model parameters
rho = 1.225; % kg.m^-3 
A = 4*3.14*0.02^2; % m^2
m = 0.00283; % kg
g = 9.81; % m.s^-2
C_D = 0.47; % -
u_eq = (2*m*g/(rho*A*C_D*alpha^2))^(1/(2*beta)); % eq. fan speed

% simulate model
load_system(model)
out = sim(model);

% extract simulation
t = out.simout(:,1);
z_sp = out.simout(:,2);
u = out.simout(:,3);
P = out.simout(:,4);
I = out.simout(:,5);
D = out.simout(:,6);
h = out.simout(:,7);
h_measured = out.simout(:,8);

% visualise simulation
figure(1)
subplot(2,1,1)
plot(t, z_sp, 'color', [0.7,0,0,1], 'LineWidth', 2, 'Displayname', 'Setpoint')
hold on
plot(t, h, 'k', 'LineStyle', '--', 'LineWidth', 2, 'Displayname', 'System (actual)')
plot(t, h_measured, 'k', 'LineWidth', 2, 'Displayname', 'System (measured)')
xline(disturb_time, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.', ...
    'DisplayName', 'Fan disturbance')
hold off
xlim([0, simtime])
%ylim([0,0.65])
ylabel('Height (m)')
title('Output')
legend('Location', 'southeast')

subplot(2,1,2)
plot(t, u, 'k', 'LineWidth', 2, 'Displayname', 'Total input')
hold on
plot(t, P, 'color', [0.7,0,0,1], 'LineWidth', 2, 'Displayname', 'P')
plot(t, I, 'color', [0,0.7,0,1], 'LineWidth', 2, 'Displayname', 'I')
plot(t, D, 'color', [0,0,0.7,1], 'LineWidth', 2, 'Displayname', 'D')
xlim([0, simtime])
ylim([-20,255])
xlabel('Time (s)')
ylabel('Fan speed (-)')
title('Input')
set(gcf,'units','inches','position',[0,0,8.3,0.40*11.7])
legend('Location', 'southeast')
% save result
saveas(gcf,'tuned_PID_insilico.pdf')