%% Open loop using Rho Prescott model
% Edited by Hubald Verzijl 
addpath('RhoPrescottModel')

clear;
close all;

% constraints on u
umin = -5;
umax = 5;

% initial state
x0 = 1;
xxx0 = x0;
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

% Optimizer settings
options = optimset('Display', 'off');

% simulation horizon
freq = 20;
N = 500*freq; % Simulation horizon (sampling freq = 20000Hz)

% Biphasic pulses stimulation 
PW = 0.25;      % 250us
PW_inter = 20;  % 50Hz
PW_inter = 2;   % 500Hz
stim = [];
for i = 1:PW_inter*freq
    if (i < PW*freq*0.5)
        stim(i,1) = -2;
    elseif (i < PW*freq)
%         stim(i,1) = 2;
    else
        stim(i,1) = 0;
    end
end
stim = repmat(stim,ceil(N/(PW_inter*freq))+1,1);

while (k <= N)
    if (k<=controller_start*freq) || (k>=controller_stop*freq)
        controller_active = false;
    else 
        controller_active = true;
    end
        
    uuukplusj = stim(k+1);
    if ~controller_active
        uuukplusj = 0*uuukplusj;
    end
    uuu0Nminus1 = [uuu0Nminus1  uuukplusj];

    % Update real model output 
    Istim = Ibase + rand(1)*20;
    xk_RP_OL = xk_RP_OL + dt*D2(i,xk_RP_OL, Istim);
    xN_RP_OL = [xN_RP_OL xk_RP_OL];
    Istim = Ibase + rand(1)*20 + uuukplusj*20;
    xk_RP_CL = xk_RP_CL + dt*D2(i,xk_RP_CL, Istim);
    xN_RP_CL = [xN_RP_CL xk_RP_CL];

    k = k + 1;                
   
    percentage = round(100*k/(N), 2);
    disp(strcat("Progress: ", num2str(percentage),'%'));
end

% Create arrays to plot
u0Nminus1 = uuu0Nminus1;

x1N_RP_OL = xN_RP_OL(1:1,:);
x1N_RP_CL = xN_RP_CL(1:1,:);

% Plot 
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
ylabel('Membrane Potential [$mV$]','Interpreter','latex','fontsize',18)

yyaxis(ax1,'right')
h(4) = plot(0:(1/freq):(length(u0Nminus1)-1)/freq,u0Nminus1,'Color',[0.07 0.53 0.03],'LineWidth',1);
h(5) = plot([0 N/freq], [umin umin],'--','Color',[0.07 0.53 0.03]);
h(6) = plot([0 N/freq], [umax umax],'--','Color',[0.07 0.53 0.03]);
ax1.YAxis(2).Color = [0 0.5 0];

leg1 = legend([h(1), h(2), h(3), h(4), h(5)],...
        'Potential without stimulation ',...
        'Potential with predetermined stimulation',...
        'Stimulus active interval',...
        'Actuation current ',...
        'Actuation current limits');
set(leg1,'Interpreter','latex','fontsize',14);
ylim([-6 6])
ylabel('Actuaion Current [$\mu A$]','Interpreter','latex','fontsize',18)

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