%% Plot for differen stimulation currents
clear; clc; close;

set(0,'defaulttextInterpreter','latex')
format = {'fontsize',18}; % name, value pairs
formatLegend = {'fontsize',14};

I_vec = [0.3,0.5,0.7];
% I = 0.5; % uA
tspan = 2500;

dt = 0.01;               % time step for forward euler method

addpath('gating')

% Time parameters

Plot = 0;
multiPlot = 1;

loop  = ceil(tspan/dt);   % no. of iterations of euler

% Ininitial values
v_init = -56.8;

% Initializing variable vectors
t   = (1:loop)*dt;

V   = zeros(loop,1);
mf  = zeros(loop,1);
hf  = zeros(loop,1);

ml1 = zeros(loop,1);
hl1 = zeros(loop,1);

ml2 = zeros(loop,1);
hl2 = zeros(loop,1);

n   = zeros(loop,1);
nnorm = zeros(loop,1);

V(1) = v_init;
mf(1) = minff(v_init);
ml1(1) = minfl1(v_init);
ml2(1) = minfl2(v_init);
hf(1) = hinff(v_init);
hl1(1) = hinfl1(v_init);
hl2(1) = 0.5;
n(1) = 0;
nnorm(1) = ninfnorm(v_init);

HH = [V(1); mf(1); hf(1); ml1(1); hl1(1); ml2(1); hl2(1); n(1); nnorm(1)];

figure('Renderer', 'painters', 'Position', [10 10 1200 300])
ax1 = nexttile;

% Euler method
for k = 1:length(I_vec)
    I = I_vec(k);
    for i=1:loop-1
        HH = HH_step(HH, dt, I);
        V(i+1) = HH(1);
    end
    
    subplot(1,length(I_vec),k)

    plot(t,V,'-','Color',[0 0.28 0.67],'LineWidth',1);
    hold on; grid on;
    set(gca,'GridLineStyle',':')
    axis([min(t) max(t) -70 10])
    xlabel('Time [ms]',format{:});
    label = append('$v_m(t)$ for $I(t)=$ ',num2str(I),' $\mu$A');
    legend(label,formatLegend{:},'interpreter','latex');
    
    if k == 1
        ylabel('Membrane Potential [mV]',format{:});
    end

end

% Plot 

HH_data = [t; V'];
%% Gating variables plot
vec = -85:60;
counter = 1;
for i = vec
    res1(counter) = minfl1(vec(counter));
    res2(counter) = hinfl1(vec(counter));
    counter = counter + 1;
end

plot(vec,res1)
hold on
plot(vec,res2)



