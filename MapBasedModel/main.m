%% Extended map based model
% Created by Hubald Verzijl 

clear; close; clc;

set(0,'defaulttextInterpreter','latex')
format = {'fontsize',18}; % name, value pairs
formatLegend = {'fontsize',14};

sigmas = [-0.005, -0.001, 0.001];
n = 1:2000;
X(:,1) = [-1; -0.1];
global alpha beta mu
alpha = 0.99;
beta = 0;
mu = 0.02;

figure('Renderer', 'painters', 'Position', [10 10 1200 300])
ax1 = nexttile;

for k = 1:length(sigmas)
    sigma = sigmas(k);
    
    for i = n
        zeta = 0.002*randn(1);
        X_bar = map_model(X(:,i), sigma, zeta);
        X(:,i+1) = X_bar;
    end
    
    output_data(k,:) = X(1,1:end-1);
    
    subplot(1,length(sigmas),k)
    plot(n, X(1,1:end-1),'-','Color',[0 0.28 0.67],'LineWidth',1);
    hold on; grid on;
    set(gca,'GridLineStyle',':')
%     plot(n, X(2,1:end-1),'Color',[1,0,0,0.5]);
    axis([min(n) max(n) -1.5 1])
    xlabel('Time [$n$]',format{:});
    label = append('$x_k$ for $\sigma=',num2str(sigmas(k)),'$');
    legend(label,formatLegend{:},'interpreter','latex');
    
    if k == 1
        ylabel('State $x_k$',format{:})
    end
end



MB_data = [n; output_data];