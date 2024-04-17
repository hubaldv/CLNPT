%% Map-based neuron
clear; close; clc;

set(0,'defaulttextInterpreter','latex')
format = {'fontsize',18}; % name, value pairs
formatLegend = {'fontsize',14};

%% Cobweb plot - Oscillations

y = -0.2;
x = -2:0.001:2;
for i = 1:length(x)
    [xi(i), d(i)] = fAlpha(x(i), y);
end

figure
% figure('Renderer', 'painters', 'Position', [10 10 1200 500])
% subplot(1,2,1)
hold on; grid on;
plot(x,x,'k-.','LineWidth',1);
plot(-1.4309,-1.4309,'r*','MarkerFaceColor','k','MarkerSize',14);

for i = 1:4
    plot(x(d==i),xi(d==i),'Color',[0 0.28 0.67],'LineWidth',2);
end

data = [-1; -1];
for i = 1:10
    x1 = data(2,end);
    xn1 = fAlpha(x1,y);
    data = [data, [x1; xn1], [xn1; xn1]];
end

plot(data(1,:), data(2,:),'k','LineWidth',1);
plot(-1,-1.1,'kv','MarkerFaceColor','k');
plot(-1.1,-1.19,'k<','MarkerFaceColor','k');
plot(-1.19,-1.26,'kv','MarkerFaceColor','k');

xlabel('$x_{k}$', format{:});
ylabel('$x_{k+1}$', format{:});
axis equal;
axis([-2, 2, -2, 2])

legend('$x_k=x_{k+1}$','Stable point','$f_\alpha(x_k,-0.2)$',formatLegend{:},'interpreter','latex')

%% Cobweb plot - Discharge

y = 0.1;
x = -2:0.001:2;
for i = 1:length(x)
    [xi(i), d(i)] = fAlpha(x(i), y);
end

figure
% subplot(1,2,2)
hold on; grid on;
plot(x,x,'k-.','LineWidth',1);

for i = 1:4
    plot(x(d==i),xi(d==i),'Color',[0 0.28 0.67],'LineWidth',2);
end

data = [-1; -1];
for i = 1:20
    x1 = data(2,end);
    xn1 = fAlpha(x1,y);
    data = [data, [x1; xn1], [xn1; xn1]];
end

plot(data(1,:), data(2,:),'k','LineWidth',1);
plot(1+y,0,'kv','MarkerFaceColor','k');
plot(0,-1,'k<','MarkerFaceColor','k');
plot((y+1)/2,1+y,'k>','MarkerFaceColor','k');
plot(.1802,0.6,'k^','MarkerFaceColor','k');

xlabel('$x_{k}$', format{:});
ylabel('$x_{k+1}$', format{:});
axis equal;
axis([-2, 2, -2, 2])

legend('$x_k=x_{k+1}$','$f_\alpha(x_k,0.2)$',formatLegend{:},'interpreter','latex')

%% Fuction
function [f, p] = fAlpha(x,y) 
    global alpha
    if (x < (-1-alpha/2))
        f = -alpha^2/4-alpha + y;
        p = 1;
    elseif (((-1-alpha/2) <= x) && (x <= 0))
        f = alpha*x+(x+1)^2+y;
        p = 2;
    elseif ((0 < x) && (x < (y+1)))
        f = y+1;
        p = 3;
    elseif x >= (y+1)
        f = -1;
        p = 4;
    else
        disp('error')
        f = 0;
        p = 5;
    end
end


