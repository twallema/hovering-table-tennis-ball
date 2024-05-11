% clear all variables/windows
clc, clear all, close all

% load data
df = readtable('dataset.csv');
d = 0.65 - df.height;
v = df.voltage;

% define the inverse of the logistic function: v = v_max - (v_max - v_min)/(1+a*exp(-b*(h-c)))
logitModel = @(v,v_min,v_max,a,b,c)(c - (1/b)*(log((v_max-v_min)./(v_max-v)-1) - log(a)));
v_min = 0.70;
c = 0;

% define a penalty function and wrap it to be compatible with fminsearch
SSE = @(v,v_min,v_max,a,b,c,d)(sum((d-logitModel(v,v_min,v_max,a,b,c)).^2));
f = @(theta)SSE(v,v_min,theta(1),theta(2),theta(3),c,d);

% initial parameters
theta_0 = [3.3, 15, 15];

% fit logit model
[theta_min,SSE_min] = fminsearch(f,theta_0);

% compute modeled values
v_hat = linspace(0.7,3.2,200);
d_hat = logitModel(v_hat,v_min,theta_min(1),theta_min(2),theta_min(3),c);

% visualise fit and save figure
f = figure('Name','Calibration IR sensor','NumberTitle','off');
plot(v_hat,62-100*d_hat,'black','LineWidth',1.5)
hold on
scatter(v,62-100*d,50,'filled','o','black')
legend('model','data','Location','southeast','FontSize',11)
hold off
xlim([0.5,3.2])
ylim([0,70])
xlabel('voltage (V)','FontSize',11)
ylabel('height (cm)','FontSize',11)
grid("off")
set(gcf,'units','inch','position',[0,0,8.3,0.25*11.7])
ax=gca;
exportgraphics(ax,'calibration_sensor.pdf')
