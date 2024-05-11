% clear all variables, close all windows
clc, clear all, close all;

% =====================
% experiment parameters
% =====================

N = 2; % number of repeated experiments at same fan speed
u_maxs = [170, 185, 255]; % uPWM's of initial step response
u_steps = [6, 10]; % step sizes around eq.
u_eq = 151; % estimated equillibrium fan speed

initial_step_length = 7.5; % s
consecutive_step_length = 2.5; % s
output_name = 'blockResponse.csv'; % name of output 

% ===================
% perform experiments
% ===================

% compute total experiment time
simtime = initial_step_length + 9*consecutive_step_length;
% pre-initialise output
data = table([], [], [], [], [], [], 'VariableNames', ...
    {'ID', 'u_max', 'u_step', 'u', 't', 'h'});
% pre-initialise dead times
t_d = [];
% perform experiments and save sampled data in structure
counter = 1;
for i = 1:length(u_maxs)
    % set fan speed
    u_max = u_maxs(i);
    % display to user
    disp(['Working on u_max = ', num2str(u_max)])
    for j = 1:length(u_steps)
        % set step size
        u_step = u_steps(j);
        % perform experiment N times
        for k = 1:N
            % run experiment
            load_system('blockResponse.slx')
            set_param('blockResponse','ConnectedIO','on');
            out=sim('blockResponse.slx');
            % extract data
            t = out.simout(:,1);
            h = out.simout(:,2);
            u = out.simout(:,3);
            % assign height to a temperorary storage matrix
            if k == 1
                h_vect = zeros(length(t),N);
                h_vect(:,1) = h;
            else
                h_vect(:,k) = h;
            end
        end
        % average measured heights
        h = mean(h_vect, 2);
        % assign to data
        new_data = table(counter*ones(length(t),1), u_max*ones(length(t),1), ...
            u_step*ones(length(t),1), u, t, h, 'VariableNames', ...
            {'ID', 'u_max', 'u_step', 'u', 't', 'h'});
        data = vertcat(data, new_data);
    % advance ID counter
    counter = counter + 1;
    end
    % show data from same u_max to user
    figure(1)
    for j = 1:length(u_steps)
        % slice data
        d = data(data.u_max==u_max & data.u_step==u_steps(j),:);
        % visualise experiment
        plot(data.t(data.u_max==u_max & data.u_step==u_steps(j)), ...
            data.h(data.u_max==u_max & data.u_step==u_steps(j)), ...
            'LineWidth', 2)
        hold on
    end
    hold off
    xlim([0,5])
    title("Approximate lift-off dead time?")
    % prompt & save dead time
    t_d(i) = input("Approximate lift-off dead time? ");
end

% =================
% visualise results
% =================

% visualise all experiments
figure(2)
for i = 1:length(u_maxs)
    % select right subplot
    subplot(length(u_maxs), 1, i); 
    for j = 1:length(u_steps)
        % plot height
        plot(data.t(data.u_max==u_maxs(i) & data.u_step==u_steps(j)), ...
            data.h(data.u_max==u_maxs(i) & data.u_step==u_steps(j)), ...
            'DisplayName',['u_{step} = ', num2str(u_steps(j))], ...
            'LineWidth', 2)
        hold on
    end
    % vertical line at dead time
    xline(t_d(i), '-.', {['t_d = ', num2str(t_d(i)),' s']}, 'Color', 'k', 'LineWidth', 1, ...
        'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'left')
    % vertical lines at steptimes
    xline(2*initial_step_length, '-.', 'Color', 'k', 'LineWidth', 1, ...
        'HandleVisibility', 'off')
    xline(2*initial_step_length+3*consecutive_step_length, '-.', 'Color', ...
        'k', 'LineWidth', 1, 'HandleVisibility', 'off')
    % set decorations
    title(['u_{max} = ', num2str(u_maxs(i))])
    ylabel('height (m)')
    xlim([0,27])
    ylim([0,0.6])
    if i == length(u_maxs)
        legend('Location', 'southeast')
    end
end
% set more decorations
xlabel('time (s)')
set(gcf,'units','inches','position',[0,0,8.3,length(u_maxs)*0.20*11.7])
% save figure
saveas(gcf,'blockResponse.pdf')

% =================
% cut out dead time 
% =================

for i = 1:length(u_maxs)
    data.t(data.u_max==u_maxs(i)) = data.t(data.u_max==u_maxs(i)) - t_d(i);
end
data = data(data.t >= 0,:);

% ============
% save results
% ============

writetable(data,output_name,'Delimiter',',')  