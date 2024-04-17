%% Create input output data
% Main file for the Rho and Prescott model
% Run from within RhoPrescottModel folder

clear; close; clc;

set(0,'defaulttextInterpreter','latex')
format = {'fontsize',18}; % name, value pairs
formatLegend = {'fontsize',14};

result = [-65; 0.3];
row = 1;

dt = 0.05;
t_max = 500;
vec = 0:dt:t_max;

stims = [31.5, 32, 32.5];

figure('Renderer', 'painters', 'Position', [10 10 1200 300])
ax1 = nexttile;

for k = 1:length(stims)
    row = 1;

    for i = 0:dt:t_max-dt
        row = row + 1;

        w = rand(1)*20;
        Istim = stims(k) + w;
    %     disp(Istim)
        result(:,row) = result(:,row-1) + dt*D2(i,result(:,row-1), Istim);
    end
    
    subplot(1,length(stims),k)

    plot(vec,result(1,:),'-','Color',[0 0.28 0.67],'LineWidth',1);
    axis([0 t_max -90 40])
    grid on; 
    set(gca,'GridLineStyle',':')
    xlabel('Time [ms]',format{:});
    label = append('$v_m(t)$ for $I(t)=',num2str(stims(k)),'\mu A$');
    legend(label,formatLegend{:},'interpreter','latex');
    
    if k == 1
        ylabel('Membrane Potential [$mV$]',format{:});
    end
end

ML_data = [vec; result(1,:)];