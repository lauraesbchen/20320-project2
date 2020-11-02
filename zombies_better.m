%% Nonquarantine
S_o = 5.13E7; % 
R_o = 0; % 
D_o = 0; % 
Z_o = 1; %could be changed depending on assumptions about how the leak impacted people  

Vb  = 0; %birth rate, 
Ksd = 0; % 
Ksz = 20/24; %per hour
Kdz = 0; %
%Kzr = 0.00008333333333;
Kzr = 0.0008333333333;

%tspan = 60*60; %minutes
tspan1 = linspace(0,100,5E1);


options = odeset('AbsTol', 1e-12, 'RelTol', 5e-12);

[t, p] = ode23s(@zombies, tspan1, [S_o, D_o, Z_o, R_o], options, [Vb, ...
    Ksd, Ksz, Kdz, Kzr]);

figure()
semilogy(t, p(:,1), '--b','LineWidth', 3)
hold on
semilogy(t, p(:,2), 'o-c', 'LineWidth', 2)
hold on
legend(["Susceptible", "Dead"])
hold off
figure()
semilogy(t, p(:,3), '--r', 'LineWidth', 3)
hold on
semilogy(t, p(:,4), 'o-k', 'LineWidth', 2)
legend(["Zombies", "Removed"])

% legend(["Susceptible", "Dead", "Zombies", "Removed"])




%% Quarantine 

S_o = 4.12E7;
R_o = 0;
D_o = 0;
Z_o = 1;
Q_o = 10E6; 

Vb  = 0;
Ksd = 0;
Ksz = 20/24;
Kdz = 0;
Kzr = 0.00008333333333;
Kqs = 0;
Ksq = 1/6;

options = odeset('AbsTol', 1e-14, 'RelTol', 1e-14);
time_span = linspace(0,100,50);
[t_q, p_q] = ode15s(@zombies_quarantine, time_span, [S_o, D_o, Z_o, R_o, Q_o], options, [Vb, ...
    Ksd, Ksz, Kdz, Kzr, Kqs, Ksq]);

figure()
semilogy(t_q, p_q(:,1), 'b', 'LineWidth', 2)
hold on
semilogy(t_q, p_q(:,2), 'c', 'LineWidth', 2)
hold on
semilogy(t_q, p_q(:,3), 'r', 'LineWidth', 2)
hold on
semilogy(t_q, p_q(:,4), 'k', 'LineWidth', 2)
hold on 
semilogy(t_q, p_q(:,5), 'o-g', 'LineWidth', 2)
legend(["Susceptible", "Dead","Zombies", "Removed", "Quarantined"])

figure()
N = p_q(:,1)+p_q(:,3)+p_q(:,4)+p_q(:,5);
plot(t_q,N)

%% Functions
% Non-Quarantine of Zombies
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
Kzr = k(5);

%assuming no birth rate
dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./(S+Z);
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./(S+Z) + Kdz.*D - Kzr.*S.*Z./(S+Z);
dydt(4,:) = Kzr.*Z.*S./(S+Z);


end

% Quarantine People from Zombies

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
Kzr = k(5);
Kqs = k(6);
Ksq = k(7);

%assuming no birth rate
dydt(1,:) = Vb - Ksd.*S - Ksz.*S.*Z./(S+Z) - Ksq.*S +Kqs.*Q;
dydt(2,:) = Ksd.*S - Kdz.*D;
dydt(3,:) = Ksz.*S.*Z./(S+Z) + Kdz.*D - Kzr.*S.*Z./(S+Z);
dydt(4,:) = Kzr.*Z.*S./(S+Z);
dydt(5,:) = Ksq.*S - Kqs.*Q;

end 



