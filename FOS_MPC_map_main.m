%% FOS-MPC using the Mab-Based model
% Edited by Hubald Verzijl 

addpath('FOS_MPC_func')
addpath('MapBasedModel')

clear;
close all;
% Setting up parameters
randn('seed',2020);

% the sample model data goes here
load('MapBasedModel/MB_data.mat');
signal = MB_data(3,:)';
t_vec = MB_data(1,:)';
fs = 1/(mean(t_vec(2:end)-t_vec(1:end-1)));

% Find autocorrelation
lags = 500;
[autocor, autocor_lag, confidence_95] = autocorrelation(signal, lags);
% autocorr(signal, 'NumLags', lags)

% Stable FOS:
start = 1;
ending = start + 1549;
% Unstable FOS:
% start = 500;
% ending = start + 1549;

figure 
tiledlayout(2,1)
nexttile;
plot(t_vec, signal,'b');
hold on; grid on;
plot([t_vec(start), t_vec(ending)], [0, 0], 'r', 'LineWidth',2)
xlabel('Time n');
ylabel('State x]');
leg = legend('Rho and Prescott model output','Interval used for FOS system ID');
set(leg,'Interpreter','latex','fontsize',12);

nexttile;
plot([0, lags],[0,0],'k');
hold on; grid on;
for i = 0:lags
    cor = autocor(i+1);
    plot([i, i], [0, cor], 'r');
    plot([i, i], [cor, cor], 'ro');
end
h1 = plot([0 lags], [confidence_95 confidence_95],'b');
plot([0 lags], -[confidence_95 confidence_95],'b');
xlabel('Lag [samples]');
ylabel('Sample autocorrelation');
leg = legend(h1, '95\% confidence interval');
set(leg,'Interpreter','latex','fontsize',12);

%% Find A and alpha
order = WT_estimator_v3(signal(start:ending)',1);
ZZZ = lin_kernel(signal(start:ending)',1,order);

A = squeeze(ZZZ{1}(:,:,1));
alpha_FOS = order;

B = 1;
Bw = 0.1;
lambda = 1;
lambdaR = 1;

nu = size(B,2);

% size of LTI truncation
p = 15;

% initial state
x0 = 1;

% sampling frequency
freq = 20;

% simulation horizon
N = 500*freq; % Simulation horizon (sampling freq = 20000Hz)

n = size(A,1);      

AAA = FOS_LTI(A,alpha_FOS,N,n,p);

BBB = B;
BBBw = Bw;
xxx0 = x0;
for k=2:p
    BBB = [BBB; zeros(n,nu)];
    BBBw = [BBBw; zeros(n,n)];
    xxx0 = [xxx0; zeros(n,1)];
end

% C = zeros(1, size(AAA,1));
% sys = ss(AAA, BBB, C, 0);

%% MPC 

% Control objective reference
ref = -1.5*ones(size(AAA,1),1);

% prediction and control horizon
P = 50;
M = 10;

[AAAA, BBBB] = predmodgen(AAA,BBB,P);

[Htotal, ~, QQ, cc] = costgen(lambda,lambdaR,AAA,BBB,AAAA,BBBB,xxx0,n,P);

% constraints on u
umin = -0.5;
umax = 0.5;

uuu0Nminus1 = [];
xxx1N_OBS = [];
xxx1N_NO_OBS = [];
xxx1N_OLD = [];
k = 0;
xxxkplusj_OBS = xxx0;
xxxkplusj_NO_OBS = xxx0;
xxxkplusj_OLD = xxx0;

% Real model 
MB_x0 = [-1; -0.1];
xN_MB_OL = []; 
xN_MB_CL = [];
xk_MB_OL = MB_x0;
xk_MB_CL = MB_x0;

dt = 0.05;          % Sample time for the model
Ibase = 33;         % Input to the model, use this parameter to add/remove activity

controller_start = 100;
controller_stop = 350;
controller_active = true; 

% Observer using Kalman filter
ccc = [1,zeros(1,p-1)];
sys = ss(AAA,BBB,ccc,0, -1);

Q = 1;
R = 1;
[kalmf,L,~] = kalman(sys,Q,R,0);

% Optimizer settings
options = optimset('Display', 'off');

global alpha beta mu
alpha = 0.99;
beta = 0;
mu = 0.02;

sigma = -0.001;

while (k <= N)
    if (k<=controller_start*freq) || (k>=controller_stop*freq)
        controller_active = false;
    else 
        controller_active = true;
    end
        
    ftotal = 2*BBBB'*AAAA'*QQ*AAAA(:,1:size(AAA,1))*AAA*(xxxkplusj_OBS-ref) + BBBB'*AAAA'*cc;
    uuu0Pminus1 = quadprog(Htotal, ftotal,[],[],[],[], umin*ones(P*nu,1), umax*ones(P*nu,1), [], options);
    for j=0:(M-1)
        uuukplusj = uuu0Pminus1((j*nu+1):((j+1)*nu));
        if ~controller_active
            uuukplusj = 0*uuukplusj;
        end
        uuu0Nminus1 = [uuu0Nminus1  uuukplusj];
        
        % Update real model output
        zeta = 0.002*randn(1);
        xk_MB_OL = map_model(xk_MB_OL, sigma, zeta);
        xN_MB_OL = [xN_MB_OL xk_MB_OL];
       
        xk_MB_CL = map_model(xk_MB_CL, sigma+uuukplusj, zeta);
        xN_MB_CL = [xN_MB_CL xk_MB_CL];
        
        % Update FOS model output, no measurements during control horizon       
        w1 = normrnd(0,1,[n,1]);
        xxxkplusj_OBS = AAA*xxxkplusj_OBS + BBB*uuukplusj + BBBw*w1;
        xxxkplusj_NO_OBS = AAA*xxxkplusj_NO_OBS + BBB*uuukplusj + BBBw*w1;
        xxxkplusj_OLD = AAA*xxxkplusj_OLD + BBBw*w1;
        
        xxx1N_OBS = [xxx1N_OBS xxxkplusj_OBS];
        xxx1N_NO_OBS = [xxx1N_NO_OBS xxxkplusj_NO_OBS];
%         if ~controller_active
%             xxx1N_OBS(:,end) = NaN;
%         end
        xxx1N_OLD = [xxx1N_OLD xxxkplusj_OLD];
        k = k + 1;                
    end
    
    for j=0:P
        uuukplusj = zeros(nu,1);
        uuu0Nminus1 = [uuu0Nminus1  uuukplusj];
        
        % Update real model output
        zeta = 0.002*randn(1);
        xk_MB_OL = map_model(xk_MB_OL, sigma, zeta);
        xN_MB_OL = [xN_MB_OL xk_MB_OL];
       
        xk_MB_CL = map_model(xk_MB_CL, sigma+uuukplusj, zeta);
        xN_MB_CL = [xN_MB_CL xk_MB_CL];
        
        % Update FOS model output, based on real output (Kalman filter)
        w1 = normrnd(0,0,[n,1]);
        xxxkplusj_OBS = AAA*xxxkplusj_OBS + L * ([1 0]*xk_MB_CL - ccc*xxxkplusj_OBS) + BBB*uuukplusj+BBBw*w1;
        xxxkplusj_NO_OBS = AAA*xxxkplusj_NO_OBS + BBB*uuukplusj + BBBw*w1;
        xxxkplusj_OLD = AAA*xxxkplusj_OLD + BBBw*w1;
        
        xxx1N_OBS = [xxx1N_OBS xxxkplusj_OBS];
        xxx1N_NO_OBS = [xxx1N_NO_OBS xxxkplusj_NO_OBS];
%         if ~controller_active
%             xxx1N_OBS(:,end) = NaN;
%         end   
        xxx1N_OLD = [xxx1N_OLD xxxkplusj_OLD];
        k = k + 1;       
    end
    
    percentage = round(100*k/(N), 2);
    disp(strcat("Progress: ", num2str(percentage),'%'));
end

% Create arrays to plot
x1N_OBS = xxx1N_OBS(1:n,:);
x1F_NO_OBS = xxx1N_NO_OBS(1:n,:);
x1N_OLD = xxx1N_OLD(1:n,:);

u0Nminus1 = uuu0Nminus1;

x1N_MB_OL = xN_MB_OL(1:n,:);
x1N_MB_CL = xN_MB_CL(1:n,:);

% Plot 
figure('units','normalized','outerposition',[0 0 1 1])
tiledlayout(2,1)

% Top plot
ax1 = nexttile;
hold on; grid on;
yyaxis(ax1,'left')
plot(0:(1/freq):(length(x1F_NO_OBS(1,:))/freq),[mean(x0) x1F_NO_OBS(1,:)],'-','Color',[1 0 0.22],'LineWidth',1);
plot(0:(1/freq):(length(x1N_OBS(1,:))/freq),[mean(x0) x1N_OBS(1,:)],'-','Color',[0 0.28 0.67 0.5],'LineWidth',1);
ylabel('FOS output [...]','fontsize',16)

yyaxis(ax1,'right')
ax1.YAxis(2).Color = [0.07 0.53 0.03];
plot(0:(1/freq):(length(u0Nminus1)-1)/freq,u0Nminus1,'Color',[0.07 0.53 0.03],'LineWidth',2);
plot([0 N/freq], [umin umin],'--','Color',[0.07 0.53 0.03])
plot([0 N/freq], [umax umax],'--','Color',[0.07 0.53 0.03])
ylim([2*umin 2*umax])
ylabel('FOS input [...]','fontsize',16)

leg1 = legend(  '$\mathrm{FOS \: output \: without \: observer \: (left \: axis)}$',...
                '$\mathrm{FOS \: output \: with \: Kalman \: observer \: (left \: axis)}$',...
                '$\mathrm{Actuation \: Signal \: (right \: axis)}$',...
                '$\mathrm{Actuation \: Signal \: Limits \: (right \: axis)}$');
set(leg1,'Interpreter','latex','fontsize',14);
str = strcat("FOS output (LTI truncation: ",num2str(p),...
    ", Prediction horizon: ", num2str(P),...
    ", Control horizon: ", num2str(M),...
    ", fs: discrete time steps)");
title(str,'fontsize',16)
xlim([0 (length(x1N_MB_OL(1,:))/freq)])
xlabel('Time n','fontsize',16)


% Bottom plot
ax2 = nexttile;
plot(0:(1/freq):(length(x1N_MB_OL(1,:))/freq),[MB_x0(1) x1N_MB_OL],'Color',[1 0 0.22]);
hold on; grid on;
plot(0:(1/freq):(length(x1N_MB_CL(1,:))/freq),[MB_x0(1) x1N_MB_CL],'Color',[0 0.28 0.67]);
xlim([0 (length(x1N_MB_OL(1,:))/freq)])

str = strcat("Model output (controller active interval: [" ,num2str(controller_start)," - ",num2str(controller_stop),"] n)");

leg1 = legend(  'Open-loop output (Real Model)',...
                'Closed-loop output (Real Model)');
set(leg1,'Interpreter','latex','fontsize',14);

title(str,'fontsize',16)
xlabel('Time n','fontsize',16)
ylabel('State x','fontsize',16)

%% Plot
figure('Renderer', 'painters', 'Position', [10 10 1200 300])
ax1 = nexttile;

yyaxis(ax1,'left')
h(3) = patch([controller_start,controller_start,controller_stop,controller_stop],[-110,40,40,-110],[0.9 1 1],'FaceAlpha',1,'EdgeColor','none');
hold on; grid on;
set(gca,'GridLineStyle',':')
h(1) = plot(0:(1/freq):(length(x1N_MB_OL(1,:))/freq),[MB_x0(1) x1N_MB_OL],'-','Color',[1 0 0.22],'LineWidth',1);
h(2) = plot(0:(1/freq):(length(x1N_MB_CL(1,:))/freq),[MB_x0(1) x1N_MB_CL],'-','Color',[0 0.28 0.67],'LineWidth',1);
xlim([0 (length(x1N_MB_OL(1,:))/freq)])
ylim([-2 1]);
% str = strcat("Model output (controller active interval: [" ,num2str(controller_start)," - ",num2str(controller_stop),"] ms)");
ax1.YAxis(1).Color = [0 0 0];

% title(str,'fontsize',16)
xlabel('Time [n]','Interpreter','latex','fontsize',18)
ylabel('Membrane potential [arb. unit]','Interpreter','latex','fontsize',18)

yyaxis(ax1,'right')
h(4) = plot(0:(1/freq):(length(u0Nminus1)-1)/freq,u0Nminus1,'Color',[0.07 0.53 0.03],'LineWidth',1);
h(5) = plot([0 N/freq], [umin umin],'--','Color',[0.07 0.53 0.03]);
h(6) = plot([0 N/freq], [umax umax],'--','Color',[0.07 0.53 0.03]);
ax1.YAxis(2).Color = [0 0.5 0];

leg1 = legend([h(1), h(2), h(3), h(4), h(5)],...
        'Open-loop potential',...
        'Closed-loop potential',...
        'Controller active interval',...
        'Actuation current ',...
        'Actuation current Limits');
set(leg1,'Interpreter','latex','fontsize',14);
ylim([-0.6 0.6])
ylabel('Actuaion current [arb. unit]','Interpreter','latex','fontsize',18)

set(gca, 'layer', 'top');
