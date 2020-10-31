%% Nonquarantine
S_o = 3.429E6;
R_o = 0;
D_o = 0;
Z_o = 1;

Vb  = ;
Ksd = ;
Ksz = ;
Kdz = ;
Kzr = ;

options = odeset('AbsTol', 1e-14, 'RelTol', 1e-14);

[t, p] = ode15s(@zombies, [0 1000], [S_o, R_o, D_o, Z_o], options, [Vb, ...
    Ksd, Ksz, Kdz, Kzr]);

figure()
plot(t, p(:,1), 'b')
hold on
plot(t, p(:,2), 'c')
hold on
plot(t, p(:,3), 'r')
hold on
plot(t, p(:,4), 'k')
legend(["Susceptible", "Dead","Zombies", "Removed"])

%% Quarantine 

S_o = 3.429E6;
R_o = 0;
D_o = 0;
Z_o = 1;
Q_o

Vb  = ;
Ksd = ;
Ksz = ;
Kdz = ;
Kzr = ;
Kqs = ;
Ksq = ;

options = odeset('AbsTol', 1e-14, 'RelTol', 1e-14);

[t_q, p_q] = ode15s(@zombies_quarantine, [0 1000], [S_o, R_o, D_o, Z_o, Q_o], options, [Vb, ...
    Ksd, Ksz, Kdz, Kzr, Kqs, Ksq]);

figure()
plot(t_q, p_q(:,1), 'b')
hold on
plot(t_q, p_q(:,2), 'c')
hold on
plot(t_q, p_q(:,3), 'r')
hold on
plot(t_q, p_q(:,4), 'k')
hold on 
plot (t_q, p_q(:,5))
legend(["Susceptible", "Dead","Zombies", "Removed", "Quarantined"])

%% Functions
function dydt = zombies(~,y,k)

% 
% System of ODEs describing a zombie outbreak
%%S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x4 [S, D, Z, R]
% k: a 1x5 [Vb, Ksd, Ksz, Kdz, Kzr] 
%                            Vb  - Birth rate
%                            Ksd - background death rate
%                            Ksz - bite rate
%                            Kdz - zombification rate
%                            Kzr - "zombie destruction" rate
%                            
% Output (output is dydt):
% A 1X4  [dS/dt, dD/dt, dZ/dt, dR/dt]
 

dydt = zeros(4,1);
S = y(1);
D = y(2);
Z = y(3);
R = y(4);

Vb = k(1);
Ksd = k(2);
Ksz = k(3);
Kdz = k(4);
Kzr = K(5);

%assuming no birth rate
dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z + Kdz.*D - Kzr.*S.*Z;
dydt(2,:) = Kzr.*Z.*S;

end 

function dydt = zombies_quarantine(~,y,k)

% 
% System of ODEs describing a zombie outbreak
%%S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x5 [S, D, Z, R, Q]
% k: a 1x7 [Kb, Ksd, Ksz, Kdz, Kzr, Kqs, Ksq] 
%                            Vb  - Birth rate
%                            Ksd - background death rate
%                            Ksz - bite rate
%                            Kdz - zombification rate
%                            Kzr - "zombie destruction" rate
%                            
% Output (output is dydt):
% A 1X5  [dS/dt, dD/dt dZ/dt, dR/dt, dQ/dt]
 

dydt = zeros(4,1);
S = y(1);
D = y(2);
Z = y(3);
R = y(4);
Q = y(5);

Vb = k(1);
Ksd = k(2);
Ksz = k(3);
Kdz = k(4);
Kzr = K(5);
Kqs = k(6);
Ksq = k(7);

%assuming no birth rate
dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z - Ksq.*S +Kqs.*Q;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z + Kdz.*D - Kzr.*S.*Z;
dydt(4,:) = Kzr.*Z.*S;
dydt(5,:) = Ksq.*S - Kqs.*Q;

end 