function MATLABsession5_Assignment
%This file contains instructions requested by MATLAB Notes - Session 4
%
%Three solutions for are provided for a decay equation inlcuding 
%1) An ananlytical solution.
%2) A numerical solution obtained using first-order explicit time-stepping.
%3) A numerical solution obtained using ODE45.

% tHalf (years) - Half life of Caesium-137
% C0 (ppb) - Initial mass fraction of Caesium-137
% lambda (1/years) - Decay coefficient
% t (years) - Time

%Define model parameters
tHalf=30;
C0=20;

%Calculate decay coefficient
lambda=log(2)/tHalf;

%Define times of interest
t=linspace(0,100,20)';

%Evaluate an analytical solution for comparison purposes
Ca=C0*exp(-lambda*t);

%Solve using first-order explicit time-stepping
[t,C]=MYsolver(@MYodefun,t,C0,lambda);

%Solve using ODE45
options=[];
[t,Code45]=ode45(@MYodefun,t,C0,options,lambda);

%Plot results
figure(1)
plot(t,Ca,t,C,'o',t,Code45,'rx')
xlabel('Time (years)')
ylabel('Caesium-137 concentration (ppb)')
legend('Analytical','First-order explicit')

%**************************************************************************

function dCdt=MYodefun(t,C,lambda)
%The first-order decay equation
dCdt=-lambda*C;

%**************************************************************************

function [t,f]=MYsolver(odefun,t,f0,p1)
% odefun - contains the name of the subfunction containing the ODE function
% t - a vector containing the times to be solved for.
% f0 - initial condition of f
% p1 - a parameter in the ODE
% f - a vector containing the solution

%Initialize solution vector
f=zeros(size(t));
%Set the initial condition
f(1)=f0;
%Step through explicit time-stepping sequence
for n=1:numel(t)-1;
    %Obtain an estimate of the derivative using the ODE function
    dfdt=odefun(t(n,1),f(n,1),p1);
    %Define time-step
    dt=t(n+1)-t(n);
    %Evaluate y for the current time-step
    f(n+1,1)=dfdt*dt+f(n,1);
end