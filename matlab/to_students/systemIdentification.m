% clear all variables, close all windows
clc, clear all, close all;

% ================================
% What is the in-flight dead time?
% ================================

t_d = 0.50; % s

% ==================
% Load & format data
% ==================

% load calibration and validation data
data_blockResponse = readtable('../solutions/blockResponse.csv');
data_wResponse = readtable('../solutions/wResponse.csv');
% place calibration data in a more easy-to-acces sequential format
datasets={};
for i = 1:length(unique(data_blockResponse.ID))
    datasets{i} = data_blockResponse(data_blockResponse.ID == i, 4:end);
end

% =================================================
% Define mathematical models and objective function
% =================================================

function y = simulate_TF1(t, u, y0, t_d, pars)

    % -------------------------------------------
    % A function to simulate a first-order system
    % -------------------------------------------
    % t: Timesteps
    % u(t): Corresponding inputs
    % y(0): Initial condition
    % pars: [K, tau, t_d]

    % define equillibrium fan speed
    u_eq = 152;
    % unpack parameters
    K = pars(1);
    tau = pars(2);
    % define input function
    u_func = @(t_interp) interp1(t,u,t_interp,"linear",u(1));
    % system equations with saturation
    function dydt = integrate(t,y,u,u_eq,K,tau,t_d)
        lower_bound = 0.0;
        upper_bound = 0.52;
        if y(1) >= upper_bound
            dydt = min(0, (1/tau)*(K*(u(t-t_d)-u_eq) - y));
        elseif y(1) <= lower_bound
            dydt = max(0, (1/tau)*(K*(u(t-t_d)-u_eq) - y));
        else
            dydt = (1/tau)*(K*(u(t-t_d)-u_eq) - y);
        end
    end
    % integrate system using ode45
    [t_int, y] = ode45(@(t,y) integrate(t,y,u_func,u_eq,K,tau,t_d),[t(1) t(end)], y0);
    % interpolate back to input timesteps
    y = interp1(t_int,y,t);
end

function y = simulate_TF2(t, u, y0, t_d, pars)

    % --------------------------------------------
    % A function to simulate a second-order system
    % --------------------------------------------
    % t: Timesteps
    % u(t): Corresponding inputs
    % y(0): Initial condition
    % pars: [K, omega_n, zeta, t_d]

    % define equillibrium fan speed
    u_eq = 152;
    % unpack parameters
    K = pars(1);
    omega_n = pars(2);
    zeta = pars(3);
    % append zero velocity to initial states
    y0 = [y0; 0];
    % define input function
    u_func = @(t_interp) interp1(t,u,t_interp,"linear",u(1));
    % system equations with saturation
    function dydt = integrate(t,y,u,u_eq,K,omega_n,zeta,t_d)
        lower_bound = 0.0;
        upper_bound = 0.52;
        if y(1) >= upper_bound
            dydt = [min(0, y(2)); -2*zeta*omega_n*y(2) - omega_n^2*y(1) + K*(u(t-t_d)-u_eq)];
        elseif y(1) <= lower_bound
            dydt = [max(0, y(2)); -2*zeta*omega_n*y(2) - omega_n^2*y(1) + K*(u(t-t_d)-u_eq)];
        else
            dydt = [y(2);  -2*zeta*omega_n*y(2) - omega_n^2*y(1) + K*(u(t-t_d)-u_eq)];
        end
    end
    % integrate system using ode45
    [t_int, y] = ode45(@(t,y) integrate(t,y,u_func,u_eq,K,omega_n,zeta,t_d),[t(1) t(end)], y0);
    % retain only position
    y = y(:,1);
    % interpolate back to input timesteps
    y = interp1(t_int,y,t);
end

function y = simulate_ODE(t, u, y0, t_d, pars)

    % -----------------------------------------------
    % A function to simulate the non-linear ODE model
    % -----------------------------------------------
    % t: Timesteps
    % u(t): Corresponding inputs
    % y(0): Initial condition
    % pars: [k, C_L, t_d]

    % hardcode system parameters
    rho = 1.225; % kg.m^-3 
    A = 4*3.14*0.02^2; % m^2
    m = 0.00283; % kg
    g = 9.81; % m.s^-2
    C_L = 0.47; % -
    % calibrated variables
    k1 = pars(1);
    k2 = pars(2);
    % append zero velocity to initial states
    y0 = [y0; 0];
    % define input function
    u_func = @(t_interp) interp1(t,u,t_interp,"linear",u(1));
    % system equations with saturation
    function dydt = integrate(t,y,u,rho,A,m,g,k1,k2,C_L,t_d)
        lower_bound = 0.0;
        upper_bound = 0.50;
        if y(1) >= upper_bound
            dydt = [min(0,y(2)); (rho*A*C_L/(2*m))*(k1*u(t-t_d)^k2-y(2))^2-g];
        elseif y(1) <= lower_bound
            dydt = [max(0,y(2)); (rho*A*C_L/(2*m))*(k1*u(t-t_d)^k2-y(2))^2-g];
        else
            dydt = [y(2); (rho*A*C_L/(2*m))*(k1*u(t-t_d)^k2-y(2))^2-g];
        end
    end
    % integrate system using ode45
    [t_int, y] = ode45(@(t,y) integrate(t,y,u_func,rho,A,m,g,k1,k2,C_L,t_d), ...
        [t(1) t(end)], y0, odeset('RelTol', 1e-9));
    % retain only position
    y = y(:,1);
    % interpolate back to input timesteps
    y = interp1(t_int,y,t);
end

function SSE = Jtheta(pars, data, model, t_d)
    
    % --------------------------------------------
    % A sum-of-squared errors for a single dataset
    % --------------------------------------------
    % pars: vector
    %   Values of the model's parameters
    % data: table
    %   Dataset, contains columns 'u' (input), 't' (time) and 'h' (height)
    % model: function
    %    Mathematical model. Must have 't', 'u' and 'y0' as first
    %    arguments, followed by the parameters

    % simulate model (always starts at 0; noisy sensor!)
    y_model = model(data.t, data.u, 0, t_d, pars);
    % compute SSE (data occasionaly contains NaN)
    SSE = sum((fillmissing(data.h,'previous') - y_model).^2);
end

function SSE_total = Jtheta_total(pars, bounds, datasets, model, t_d)

    % --------------------------------------------------------------
    % A sum-of-squared errors function summing over all datasets and
    % restricting parameters within certain bounds
    %---------------------------------------------------------------

    SSE_total = 0;
    % restrict parameters within bounds
    for i = 1:length(bounds)
        if pars(i) <= bounds{i}(1)
            SSE_total = SSE_total + 1e9;
            pars(i) = bounds{i}(1);
        elseif pars(i) >= bounds{i}(2)
            SSE_total = SSE_total + 1e9;
            pars(i) = bounds{i}(2);
        end
    end
    % compute and sum SSEs of all datasets
    for i = 1:length(datasets)
        SSE_total = SSE_total + Jtheta(pars, datasets{i}, model, t_d);
    end
end

% ===================
% Estimate parameters
% ===================

par0 = [0.268, 0.557];
bounds = {[0, 1], [0, 2]};
par_ODE = fminsearch(@Jtheta_total,par0,[optimset('PlotFcns',@optimplotfval)],bounds,datasets,@simulate_ODE,t_d);
% print equillibrium fan speed
u_bar = (2*0.00283*9.81/(1.225*4*3.14*0.02^2*0.47*par_ODE(1)^2))^(1/(2*par_ODE(2)));

% =================
% Visualise results
% =================

u_maxs = unique(data_blockResponse.u_max);
u_steps = unique(data_blockResponse.u_step);

figure(1)
% block response
for i = 1:length(u_maxs)
    % choose correct subplot
    subplot(length(u_maxs)+1, 1, i);
    for j = 1:length(u_steps)
        % slice dataset
        data = data_blockResponse(data_blockResponse.u_max == u_maxs(i) & ...
        data_blockResponse.u_step == u_steps(j),:);
        % Plot entire timeseries
        % data
        plot(data.t, data.h, 'color', [0, 0, 0, 0.3], 'LineWidth', ...
            2, 'DisplayName', ['u_{step} = ', num2str(u_steps(j))])
        hold on
        % ODE model
        y_model = simulate_ODE(data.t, data.u, 0, t_d, par_ODE);
        plot(data.t, y_model, '--', 'Color', 'r', 'LineWidth', 2, ...
            'DisplayName','ODE model')
    end
    hold off
    % decorations
    xlim([0,27])
    ylim([0,0.60])
    ylabel('height (m)')
    title(['Calibration (u_{max} = ',num2str(u_maxs(i)),')'])
end
% validation
subplot(length(u_maxs)+1,1,length(u_maxs)+1)
% plot data + compute MAE
MAE_ODE = [];
for i = 1:length(unique(data_wResponse.ID))
    % slice dataset
    data = data_wResponse(data_wResponse.ID == i,:);
    % plot height data
    if i == 1
        plot(data.t, data.h, 'color', [0, 0, 0, 0.2], 'LineWidth', 2, ...
            'DisplayName', 'Data')
    else
        plot(data.t, data.h, 'color', [0, 0, 0, 0.2], 'LineWidth', 2, ...
            'HandleVisibility','off')
    end
    hold on
    % ODE model
    y_model = simulate_ODE(data.t, data.u, 0, t_d, par_ODE);
    MAE_ODE(i) = 100*mean(abs(fillmissing(data.h, 'previous') - y_model));
end
% plot model predictions
% ODE model
y_model = simulate_ODE(data.t, data.u, 0, t_d, par_ODE);
plot(data.t, y_model, '--', 'Color',  'r', 'LineWidth', 2, 'DisplayName', ...
    ['Model (MAE = ',num2str(mean(MAE_ODE),'%.1f'),' cm)'])    
hold off
% decorations
xlim([0,27])
ylim([0,0.60])
xlabel('time (s)')
ylabel('height (m)')
set(gcf,'units','inches','position',[0,0,8.3,0.15*(length(u_maxs)+1)*11.7])
legend('Location','southeast')
title('Validation')
% save result
saveas(gcf,'calibration.pdf')
