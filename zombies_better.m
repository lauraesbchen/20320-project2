
S_o = 3.429E6;
R_o = 0;
D_o = 0;
Z_o = 1;

Alpha = 0.005;
Beta = 0.01;
Delta = 0.0001;
Gamma = 0.01;

options = odeset('AbsTol', 1e-14, 'RelTol', 1e-14);

[t, p] = ode15s(@zombies, [0 1000], [S_o, R_o, D_o, Z_o], options, [Alpha, ...
    Beta, Delta, Gamma]);

figure()
plot(t, p(:,1), 'b')
hold on
plot(t, p(:,2), 'c')
hold on
plot(t, p(:,3), 'r')
hold on
plot(t, p(:,4), 'k')
legend(["Susceptible", "Dead","Zombies", "Removed"])

function dydt = zombies(~,y,k)

% 
% System of ODEs describing a zombie outbreak
%%S is susceptible population, Z is zombie population, R is removed population
% Inputs:
% y: a 1x4 [S, D, Z, R]
% k: a 1x4 [Alpha,Beta,Delta, Gamma] 
%                            a - alpha value in model: "zombie destruction" rate
%                            b - beta value in model: bite rate
%                            d - delta value in model: background death rate
%                            g - zombification rate
% Output (output is dydt):
% A 1X4  [dS/dt, dD/dt dZ/dt, dR/dt]
 

dydt = zeros(4,1);
S = y(1);
D = y(2);
R = y(3);
Z = y(4);

a = k(1);
b = k(2);
d = k(3);
g = k(4);

%assuming no birth rate
dydt(1,:) = -d.*S -b.*S;
dydt(2,:) = d.*S - g.*D;
dydt(3,:) = -a.*Z + b.*S + g.*D;
dydt(2,:) = a.*Z;

end 