%% Plot hl2 (inactivation h late 2)

clear; clc; close;

set(0,'defaulttextInterpreter','latex')
format = {'fontsize',18}; % name, value pairs
formatLegend = {'fontsize',14};

addpath('gating');

I = 0.5;  % uA
tspan = 5000;

dt = 0.01;               % time step for forward euler method

% Ininitial values
v = -56.8;
mfi = minff(v);
ml1i = minfl1(v);
ml2i = minfl2(v);
hfi = hinff(v);
hl1i = hinfl1(v);
hl2i = 0.5;
ni = 0;
nnormi = ninfnorm(v);

Plot = 0;

loop  = ceil(tspan/dt);   % no. of iterations of euler

% Fast Na _Shh_6
% Na fast
gnaf = 25;
% Leak
gl = 1.42;  % Leak conductance


% Intermediate Na _inter
% Na late1
gnal1 = 27;

% Slow Na _ls
% Na late2
gnal2 = 1.28e-1;

% Kalium
gk = 0*0.0026;

% Voltages (Voltages + Vrest)
Ena = 62;
Ek = -94;
El = -65.5; 

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
nnorm   = zeros(loop,1);

% Set initial values for the variables

V(1)   = v;
mf(1)  = mfi;
hf(1)  = hfi;

ml1(1) = ml1i;
hl1(1) = hl1i;

ml2(1) = ml2i;
hl2(1) = hl2i;

n(1)   = ni;
nnorm(1) = nnormi;




% Euler method
for i=1:loop-1
    % Fast Na
    Inaf = gnaf * mf(i)^3*hf(i)*(V(i)-Ena);

    % Inter Na    
    Inal1 = gnal1 * ml1(i)*hl1(i)*(V(i)-Ena);

    % Slow Na
    Inal2 = gnal2 * ml2(i)*hl2(i)*(V(i)-Ena);

    % Potassium 
    Ik = gk * n(i)*(V(i)-Ek);
    Iknorm = gk * nnorm(i)^4*(V(i)-Ek);

    % Leak
    Il = gl * (V(i)-El);

    % Update membrane potential
    V(i+1) = V(i) + dt*(-Inaf-Iknorm-Il-Inal1-Inal2-Ik+I);

    % Update ion channels dynamics
    % Fast Na
    mf(i+1) = mf(i) + dt*((minff(V(i))-mf(i))/taumf(V(i)));
    hf(i+1) = hf(i) + dt*((hinff(V(i))-hf(i))/tauhf(V(i)));

    % Inter Na
    ml1(i+1) = minfl1(V(i));
    hl1(i+1) = hl1(i) + dt*((hinfl1(V(i))-hl1(i))/tauhl1(V(i)));

    % Slow Na
    ml2(i+1) = minfl2(V(i));
    hl2(i+1) = hl2(i) + dt*((hinfl2(V(i))-hl2(i))/tauhl2(V(i)));

    % K
    n(i+1) = n(i) + dt*((ninfk(V(i))-n(i))/taunk(V(i)));    
    nnorm(i+1) = nnorm(i) + dt*((ninfnorm(V(i))-nnorm(i))/taunnorm(V(i)));
end
  
figure('Renderer', 'painters', 'Position', [10 10 1200 300])

subplot(2,1,1);
plot(t,V,'Color',[0 0.28 0.67]);
hold on; grid on;
axis([min(t) max(t) -70 10])
xlabel('Time [ms]',format{:});
label = append('$v_m(t)$ for $I(t)= ',num2str(I),' \mu$A');
legend(label,formatLegend{:},'interpreter','latex');

ylabel('Voltage [mV]',format{:});

subplot(2,1,2);
plot(t,hl2,'Color',[1 0 0.22]);
grid on;
axis([min(t) max(t) 0.4 0.6]);
ylabel('$h_{S}(t)$',format{:},'Color',[1 0 0.22]);
xlabel('Time [ms]',format{:});
    

if Plot == 1
%     figure('Renderer', 'painters', 'Position', [10 10 1200 300])
    subplot(4,2,1);
    plot(t,V);
    title('Voltage');
    subplot(4,2,3);
    plot(t,nnorm);
    title('Fast n');
    subplot(4,2,5);
    plot(t,mf);
    title('Fast m');
    subplot(4,2,7);
    plot(t,hf);
    title('Fast h');
    
    subplot(4,2,2);
    plot(t,ml1);
    title('Inter m');
    subplot(4,2,4);
    plot(t,hl1);
    title('inter h');
    subplot(4,2,6);
    plot(t,ml2);
    title('Slow m');
    subplot(4,2,8);
    plot(t,hl2);
    title('Slow h');
%     xlabel('Time');
%     ylabel('Membrane Potential');
%     title('Voltage time series');
end
