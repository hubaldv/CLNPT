function dy = D2(t,y,Istim)
% File to update the Rho and Prescott Model

V = y(1);
w = y(2);

betaw = -13; % Adjustable?
betam = -1.2;
gammaw = 10;
gammam = 18;
gfast = 20;
gslow = 20;
gleak = 2;
Ena = 50;
Ek = -100;
Eleak = -70;
C = 2;
psiw = 0.15;

minf = 0.5*(1+tanh((V-betam)/gammam));
winf = 0.5*(1+tanh((V-betaw)/gammaw));
tauw = 1/(cosh((V-betaw)/(2*gammaw)));

n = rand(1)*20*0;

dV = (-gfast*minf*(V-Ena)-gslow*w*(V-Ek)-gleak*(V-Eleak)+Istim+n)/C;
dw = psiw*(winf-w)/tauw;

dy = [dV; dw];

end

