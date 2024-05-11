% clear all variables, close all windows
clc, clear all, close all;

% choose between manual or LQR pole placement
pole_placement = 'LQR'; % alternative: "LQR"

% Define Q (LQR pole placement)
Q = [14000 0; 0 2000];

% Define K and L (manual pole placement)
K = [118, 34]; % state feedback gain matrix
L = [100; 250000]; % state observer corrector matrix

% calibrated system parameters
t_d = 0.5; % s; dead time
alpha = 0.1885;
beta = 0.6272;

% initialise experimental parameters
simtime = 15; % experiment time
z_0 = [0, 0]; % initial position and speed
z_sp = 0.4; % desired setpoint (height)
steptime = 5; % time of setpoint change
disturb_time = 20; % time of fan effectiveness manipulation
disturb_magnitude = -7.5; % percentage change in parameter alpha per unit time

% noise and filter parameters
tau = 0.05; % s; measurement RC filter time constant
noise_power = 0; %0.0008; 
noise_sample_time = 0.05;

% select model to use
model = 'ODEModelLQR.slx';

% ODE model parameters
rho = 1.225; % kg.m^-3 
A = 4*3.14*0.02^2; % m^2
m = 0.00283; % kg
g = 9.81; % m.s^-2
C_D = 0.47; % -
u_bar = (2*m*g/(rho*A*C_D*alpha^2))^(1/(2*beta)); % eq. fan speed

% define linearised system
u = u_bar;
v = 0;
a = [0, 1; 0, (rho*A*C_D/m)*(v - alpha*u^beta)];
b = [0; alpha*beta*rho*A*C_D/m*u^(beta - 1)*(alpha*u^beta - v)];
c = [1, 0];
d = [0];

% compute optimal gain matrix
if strcmp(pole_placement, 'LQR')
    K = lqr(a,b,Q,1);
end

% simulate model
load_system(model);
out = sim(model);

% extract simulation
t = out.simout(:,1);
z_sp = out.simout(:,2);
u = out.simout(:, 4);
z = out.simout(:,5);
z_measured = out.simout(:,6);
z_hat = out.simout(:,7);
v_hat = out.simout(:,8);

% visualise output and input
figure(1)
subplot(2,1,1)
plot(t, z_sp, 'color', [0.7,0,0,1], 'LineWidth', 2, 'Displayname', 'Setpoint')
hold on
plot(t, z, 'k', 'LineStyle', '--', 'LineWidth', 2, 'Displayname', 'System (actual)')
plot(t, z_measured, 'k', 'LineWidth', 2, 'Displayname', 'System (measured)')
%xline(disturb_time, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-.', ...
%    'DisplayName', 'Fan disturbance')
hold off
xlim([0, simtime])
ylim([0,0.65])
ylabel('Height (m)')
title('Output')
legend('Location', 'southeast')

subplot(2,1,2)
plot(t, u, 'k', 'LineWidth', 2, 'Displayname', 'Total input')
xlim([0,simtime])
ylim([0,255])
xlabel('Time (s)')
ylabel('Fan speed (-)')
title('Input')

set(gcf,'units','inches','position',[0,0,8.3,0.40*11.7])
legend('Location', 'southeast')
% save result
saveas(gcf,'tuned_LQR_insilico.pdf')

% visualise state observer estimates
figure(2)
subplot(2,1,1)
plot(t, z_measured, 'k', 'LineWidth', 1, 'Displayname', 'Filtered measurement')
hold on
plot(t, z_hat, 'r', 'LineStyle', '--', 'LineWidth', 2, 'Displayname', 'Observer estimate')
hold off
ylim([0,0.65])
ylabel('Height (m)')
yyaxis right
plot(t, z_measured-z_hat, 'Color', [1,0.7,0,1], 'LineWidth', 1, 'Displayname', 'Observer error')
ylabel('Observer error (m)')
xlim([0, simtime])
title('Height')
legend('Location', 'southeast')

subplot(2,1,2)
plot(t, v_hat, 'r', 'LineStyle', '--', 'LineWidth', 2, 'Displayname', 'Observer estimate')
hold on
yline(0, 'k', 'LineWidth', 1, 'Displayname', 'v = 0 m/s')
hold off
ylabel('Velocity (m/s)')
xlim([0, simtime])
title('Velocity')
legend('Location', 'southeast')

set(gcf,'units','inches','position',[0,0,8.3,0.40*11.7])
legend('Location', 'northeast')
% save result
saveas(gcf,'state_observer.pdf')
