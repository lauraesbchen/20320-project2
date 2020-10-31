zombies1(0.0005, 10, 0.0001, 100, 1)

function [ ] = zombies1(a,b,d,T,dt)
    % This function will solve the system of ODE’s for the basic model used in
    % the Zombie Dynamics project for MAT 5187. It will then plot the curve of
    % the zombie population based on time.
    % Function Inputs: a - alpha value in model: "zombie destruction" rate
    %                  b - beta value in model: "new zombie" rate
    %                  ze - zeta value in model: zombie resurrection rate
    %                  No zombie resurrection in Train to Butan
    %              
    %                  d - delta value in model: background death rate
    %                  T - Stopping time
    %                  dt - time step for numerical solutions
    % Created by Philip Munz, November 12, 2008
    %Initial set up of solution vectors and an initial condition

    N = 5.164E7; %N is the population
    n = T/dt;
    t = zeros(1,n+1);
    s = zeros(1,n+1);
    z = zeros(1,n+1);
    r = zeros(1,n+1);
    s(1) = N;
    z(1) = 1;
    r(1) = 0;
    t = 0:dt:T;
    % Define the ODE’s of the model and solve numerically by Euler’s method:
    for i = 1:n
        s(i+1) = s(i) + dt*(-b*s(i)*z(i)); %here we assume birth rate
    %= background deathrate, so only term is -b term
        z(i+1) = z(i) + dt*(b*s(i)*z(i) -a*s(i)*z(i));
        r(i+1) = r(i) + dt*(a*s(i)*z(i) +d*s(i));
        if s(i)<0 || s(i) >N
            break
        end
        if z(i) > N || z(i) < 0
            break
        end
        if r(i) <0 || r(i) >N
            break
        end
    end
    
    hold on
    plot(t,s,'b');
    plot(t,z,'r');
    legend('Suscepties','Zombies')
end



function [z] = eradode(a,b,ze,d,Ti,dt,s1,z1,r1)
    % This function will take as inputs, the initial value of the 3 classes.
    % It will then apply Eulers method to the problem and churn out a vector of
    % solutions over a predetermined period of time (the other input).
    % Function Inputs: s1, z1, r1 - initial value of each ODE, either the
    % actual initial value or the value after the
    % impulse.
    % Ti - Amount of time between inpulses and dt is time step
    % Created by Philip Munz, November 21, 2008
    k = Ti/dt;
    %s = zeros(1,n+1);
    %z = zeros(1,n+1);
    %r = zeros(1,n+1);
    %t = 0:dt:Ti;
    s(1) = s1;
    z(1) = z1;
    r(1) = r1;
    for i=1:k
        s(i+1) = s(i) + dt*(-b*s(i)*z(i)); %here we assume birth rate
        %= background deathrate, so only term is -b term
        z(i+1) = z(i) + dt*(b*s(i)*z(i) -a*s(i)*z(i) +ze*r(i));
        r(i+1) = r(i) + dt*(a*s(i)*z(i) +d*s(i) - ze*r(i));
    end
    plot(t,z)
end

function [] = erad(a,b,ze,d,k,T,dt)
    % This is the main function in our numerical impulse analysis, used in
    % conjunction with eradode.m, which will simulate the eradication of
    % zombies. The impulses represent a coordinated attack against zombiekind
    % at specified times.
    % Function Inputs: a - alpha value in model: "zombie destruction" rate
    %                  b - beta value in model: "new zombie" rate
    %                  ze - zeta value in model: zombie resurrection rate
    %                  d - delta value in model: background death rate
    %                  k - "kill" rate, used in the impulse
    %                  T - Stopping time
    %                  dt - time step for numerical solutions
    % Created by Philip Munz, November 21, 2008

    N = 1000;
    Ti = T/4; %We plan to break the solution into 4 parts with 4 impulses
    n = Ti/dt;
    m = T/dt;
    s = zeros(1,n+1);
    z = zeros(1,n+1);
    r = zeros(1,n+1);
    sol = zeros(1,m+1); %The solution vector for all zombie impulses and such
    t = zeros(1,m+1);
    s1 = N;
    z1 = 0;
    r1 = 0;
    %i=0; %i is the intensity factor for the current impulse
    %for j=1:n:T/dt
    % i = i+1;
    % t(j:j+n) = Ti*(i-1):dt:i*Ti;
    % sol(j:j+n) = eradode(a,b,ze,d,Ti,dt,s1,z1,r1);
    % sol(j+n) = sol(j+n)-i*k*sol(j+n);
    % s1 = N-sol(j+n);
    % z1 = sol(j+n+1);
    % r1 = 0;
    %end
    sol1 = eradode(a,b,ze,d,Ti,dt,s1,z1,r1);
    sol1(n+1) = sol1(n+1)-1*k*sol1(n+1); %347.7975;
    s1 = N-sol1(n+1);
    z1 = sol1(n+1);
    r1 = 0;
    sol2 = eradode(a,b,ze,d,Ti,dt,s1,z1,r1);
    sol2(n+1) = sol2(n+1)-2*k*sol2(n+1);
    s1 = N-sol2(n+1);
    z1 = sol2(n+1);
    r1 = 0;
    sol3 = eradode(a,b,ze,d,Ti,dt,s1,z1,r1);
    sol3(n+1) = sol3(n+1)-3*k*sol3(n+1);
    s1 = N-sol3(n+1);
    z1 = sol3(n+1);
    r1 = 0;
    sol4 = eradode(a,b,ze,d,Ti,dt,s1,z1,r1);
    sol4(n+1) = sol4(n+1)-4*k*sol4(n+1);
    s1 = N-sol4(n+1);
    z1 = sol4(n+1);
    r1 = 0;
    sol=[sol1(1:n),sol2(1:n),sol3(1:n),sol4];
    t = 0:dt:T;
    t1 = 0:dt:Ti;
    t2 = Ti:dt:2*Ti;
    t3 = 2*Ti:dt:3*Ti;
    t4 = 3*Ti:dt:4*Ti;
    %plot(t,sol)
    hold on
    plot(t1(1:n),sol1(1:n),'k')
    plot(t2(1:n),sol2(1:n),'k')
    plot(t3(1:n),sol3(1:n),'k')
    plot(t4,sol4,'k')
    hold off
end