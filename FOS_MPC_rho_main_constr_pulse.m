%% Obtain the FOS using raw data
% Start with data input
% Edited by Hubald Verzijl 
addpath('FOS_MPC_func')
addpath('RhoPrescottModel')

clear;
close all;
% Setting up parameters
randn('seed',2020);

% the sample model data goes here
load('RhoPrescottModel/ML_data.mat');
signal = ML_data(2,5000:20000)';
t_vec = ML_data(1,5000:20000)'-ML_data(1,5000);
fs = 1/(mean(t_vec(2:end)-t_vec(1:end-1))/1000);

% Find autocorrelation
lags = 100;
[autocor, autocor_lag, confidence_95] = autocorrelation(signal, lags);
% autocorr(signal, 'NumLags', lags)

start = 7500;
ending = start + 549;
% start = 500;
% ending = start + 1549;

figure 
tiledlayout(2,1)
nexttile;
plot(t_vec, signal,'b');
hold on; grid on;
plot([t_vec(start), t_vec(ending)], [0, 0], 'r', 'LineWidth',2)
xlabel('Time [ms]');
ylabel('Potential [mV]');
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
alpha = order;

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

AAA = FOS_LTI(A,alpha,N,n,p);

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
ref = -40*ones(size(AAA,1),1);

% prediction and control horizon
P = 50;
M = 20;

[AAAA, BBBB] = predmodgen(AAA,BBB,P);

[Htotal, ~, QQ, cc] = costgen(lambda,lambdaR,AAA,BBB,AAAA,BBBB,xxx0,n,P);

% constraints on u
umin = -5;
umax = 5;

uuu0Nminus1 = [];
xxx1N_OBS = [];
xxx1N_NO_OBS = [];
xxx1N_OLD = [];
k = 0;
xxxkplusj_OBS = xxx0;
xxxkplusj_NO_OBS = xxx0;
xxxkplusj_OLD = xxx0;

% Real model 
RP_x0 = [-65; 0.3];
xN_RP_OL = []; 
xN_RP_CL = [];
xk_RP_OL = RP_x0;
xk_RP_CL = RP_x0;

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

uuu0Options = [];
uuu0fval = [];

while (k <= N)
    if (k<=controller_start*freq) || (k>=controller_stop*freq)
        controller_active = false;
    else 
        controller_active = true;
        
        ftotal = 2*BBBB'*AAAA'*QQ*AAAA(:,1:size(AAA,1))*AAA*(xxxkplusj_OBS-ref) + BBBB'*AAAA'*cc;
%         uuu0Pminus1 = quadprog(Htotal, ftotal,[],[],[],[], umin*ones(P*nu,1), umax*ones(P*nu,1), [], options);

        % Nonlinear optimization
        fun = @(x) 1/2*x'*Htotal*x + ftotal'*x;
        for i = 1:floor(M/2)
            [uuu0Options(:,i),uuu0fval(i)] = fmincon(fun, zeros(P,1),[],[],[],[], umin*ones(P*nu,1), umax*ones(P*nu,1), @(x)constr_pulse(x,i), options);
        end
        [~,I] = min(uuu0fval);
        uuu0Pminus1 = uuu0Options(:,I);
    end
        
    for j=0:(M-1)
        if ~controller_active
            uuukplusj = 0;
        else 
            uuukplusj = uuu0Pminus1((j*nu+1):((j+1)*nu));
        end
        uuu0Nminus1 = [uuu0Nminus1  uuukplusj];
        
        % Update real model output
        w2 = rand(1)*20;  
        Istim = Ibase + w2;
        xk_RP_OL = xk_RP_OL + dt*D2(i,xk_RP_OL, Istim);
        xN_RP_OL = [xN_RP_OL xk_RP_OL];
        Istim = Ibase + uuukplusj*20 + w2;
        xk_RP_CL = xk_RP_CL + dt*D2(i,xk_RP_CL, Istim);
        xN_RP_CL = [xN_RP_CL xk_RP_CL];
        
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
        w2 = rand(1)*20;
        Istim = Ibase + w2;
        xk_RP_OL = xk_RP_OL + dt*D2(i,xk_RP_OL, Istim);
        xN_RP_OL = [xN_RP_OL xk_RP_OL];
        Istim = Ibase + uuukplusj*20 + w2;
        xk_RP_CL = xk_RP_CL + dt*D2(i,xk_RP_CL, Istim);
        xN_RP_CL = [xN_RP_CL xk_RP_CL];
        
        % Update FOS model output, based on real output (Kalman filter)
        w1 = normrnd(0,1,[n,1]);
        xxxkplusj_OBS = AAA*xxxkplusj_OBS + L * ([1 0]*xk_RP_CL - ccc*xxxkplusj_OBS) + BBB*uuukplusj+BBBw*w1;
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

x1N_RP_OL = xN_RP_OL(1:n,:);
x1N_RP_CL = xN_RP_CL(1:n,:);

%% Plot 
% Bottom plot
fig = figure('Renderer', 'painters', 'Position', [10 10 1200 300]);
ax1 = nexttile;

yyaxis(ax1,'left')
h(3) = patch([controller_start,controller_start,controller_stop,controller_stop],[-110,40,40,-110],[0.9 1 1],'FaceAlpha',1,'EdgeColor','none');
hold on; grid on;
set(gca,'GridLineStyle',':')
h(1) = plot(0:(1/freq):(length(x1N_RP_OL(1,:))/freq),[RP_x0(1) x1N_RP_OL],'-','Color',[1 0 0.22],'LineWidth',1);
h(2) = plot(0:(1/freq):(length(x1N_RP_CL(1,:))/freq),[RP_x0(1) x1N_RP_CL],'-','Color',[0 0.28 0.67],'LineWidth',1);
xlim([0 (length(x1N_RP_OL(1,:))/freq)])
ylim([-80 40]);
% str = strcat("Model output (controller active interval: [" ,num2str(controller_start)," - ",num2str(controller_stop),"] ms)");
ax1.YAxis(1).Color = [0 0 0];

% title(str,'fontsize',16)
xlabel('Time [ms]','Interpreter','latex','fontsize',18)
ylabel('Membrane potential [mV]','Interpreter','latex','fontsize',18)

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
        'Actuation current limits');
set(leg1,'Interpreter','latex','fontsize',14);
ylim([-6 6])
ylabel('Actuaion current [$\mu$A]','Interpreter','latex','fontsize',18)

set(gca, 'layer', 'top');

% %% Bode plot
% % signal = double(y);
% 
% L = length(t_vec) ;
% Y = fft(signal);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% % P1 = 20*log(P1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% 
% semilogx(f,P1, 'Color', [0, 0, 1, 1]);
% grid on;
% title('Single-Sided Amplitude Spectrum of Y(t)')
% xlabel('f (Hz)')
% ylabel('|Y(jw)| [uV]')
% % ylabel('dB')
% % axis([0 1000 -inf inf])

