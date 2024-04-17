function [HH] = HH_step(HH, dt, I)
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

    V = HH(1); 
    mf = HH(2);
    hf = HH(3); 
    ml1 = HH(4); 
    hl1 = HH(5); 
    ml2 = HH(6); 
    hl2 = HH(7); 
    n = HH(8); 
    nnorm = HH(9);
    
    % Fast Na
    Inaf = gnaf * mf^3*hf*(V-Ena);

    % Inter Na    
    Inal1 = gnal1 * ml1*hl1*(V-Ena);

    % Slow Na
    Inal2 = gnal2 * ml2*hl2*(V-Ena);

    % Potassium 
    Ik = gk * n*(V-Ek);
    Iknorm = gk * nnorm^4*(V-Ek);

    % Leak
    Il = gl * (V-El);
    
    % Update membrane potential
    dV = -Inaf-Iknorm-Il-Inal1-Inal2-Ik+I;
    V = V + dt*dV;
    
    % Update ion channels dynamics
    % Fast Na
    mf = mf + dt*((minff(V)-mf)/taumf(V));
    hf = hf + dt*((hinff(V)-hf)/tauhf(V));

    % Inter Na
    ml1 = minfl1(V);
    hl1 = hl1 + dt*((hinfl1(V)-hl1)/tauhl1(V));

    % Slow Na
    ml2 = minfl2(V);
    hl2 = hl2 + dt*((hinfl2(V)-hl2)/tauhl2(V));

    % K
    n = n + dt*((ninfk(V)-n)/taunk(V));    
    nnorm = nnorm + dt*((ninfnorm(V)-nnorm)/taunnorm(V));
    
    HH = [V; mf; hf; ml1; hl1; ml2; hl2; n; nnorm];
end



