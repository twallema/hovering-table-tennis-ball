% clear all variables, close all windows
clc, clear all, close all;

% =====================
% experiment parameters
% =====================

N = 6; % number of runs
initial_step_length = 7.5; % s
consecutive_step_length = 2.5; % s
u_max = 170; % input u of initial step
u_eq = 152.5; % estimated equillibrium input u
u_step = 8; % step size (units u above or below u_eq)
output_name = 'wResponse.csv'; % name of output

% ===================
% perform experiments
% ===================

% compute total experiment time
simtime = 2*initial_step_length + 6*consecutive_step_length;
% pre-initialise output
data = table([], [], [], [], 'VariableNames', {'ID', 'u', 't', 'h'});

for i = 1:N
    % run experiment
    load_system('wResponse.slx')
    set_param('wResponse','ConnectedIO','on');
    out=sim('wResponse.slx');
    % assign to results
    new_data = table(i*ones(length(out.simout(:,1)),1), out.simout(:,3), out.simout(:,1), ...
        out.simout(:,2), 'VariableNames', {'ID', 'u', 't', 'h'});
    data = vertcat(data, new_data);
end

% ===================
% determine dead time
% ===================

% visualise
figure(1)
for i = 1:N
    plot(data.t(data.ID==i), data.h(data.ID==i), 'LineWidth', 2)
    hold on
end
hold off 
xlim([0,5])
title("Approximate lift-off dead time?")
% prompt
t_d = input("Approximate lift-off dead time? ");
% cut off
data.t = data.t - t_d;
data = data(data.t >= 0,:);

% ==========================
% visualise and save results
% ==========================

% visualise data
figure(1);
for i = 1:N
    plot(data.t(data.ID==i), data.h(data.ID==i), 'Color', [0,0,0,0.3], ...
        'LineWidth', 2)
    hold on
end
ylabel('height (m)')
yyaxis right
plot(data.t((data.ID==1)), data.u(data.ID==1), 'r', 'LineWidth', 2)
ylabel('input (-)')
xlabel('time (s)')
title(['u_{max} = ', num2str(u_max)])
set(gcf,'units','inches','position',[0,0,8.3,0.20*11.7])
% all axes black
ax=gca();
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
% save figure
saveas(gcf,'wResponse.pdf')

% save result
writetable(data,output_name,'Delimiter',',')  