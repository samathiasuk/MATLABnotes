function MATLABsession4_Assignment
%This file contains instructions requested by MATLAB Notes - Session 4

%Define range of interest
x=linspace(0,2,10);

%Define constraint for integration
x0=0;
F0=1;

figure(1)
clf

%Question 1
subplot(1,3,1)
y=5*x.^2+2*x+2;
dydxAna=10*x+2;
FAna=5/3*x.^3+x.^2+2*x+1;
%Do approximate calculus
[dydx,xMP,F]=ApproxCalculus(x,y,x0,F0);
%Plot results
PlotResults(x,y,dydxAna,FAna,xMP,dydx,F)
title('1) y = 5x^2 + 2x + 2')

%Question 2
subplot(1,3,2)
y=exp(x);
dydxAna=exp(x);
FAna=exp(x);
%Do approximate calculus
[dydx,xMP,F]=ApproxCalculus(x,y,x0,F0);
%Plot results
PlotResults(x,y,dydxAna,FAna,xMP,dydx,F)
title('2) y = e^x')

%Question 3
subplot(1,3,3)
y=1./(1+x.^2);
dydxAna=-2*x./(1+x.^2).^2;
FAna=atan(x)+1;
%Do approximate calculus
[dydx,xMP,F]=ApproxCalculus(x,y,x0,F0);
%Plot results
PlotResults(x,y,dydxAna,FAna,xMP,dydx,F)
title('3) y = 1/(1+x^2)')

function PlotResults(x,y,dydxAna,FAna,xMP,dydx,F)
plot(x,y,'g',xMP,dydx,'ob',x,F,'sr',x,dydxAna,'b',x,FAna,'r')
legend('y','dy/dx finite difference','\int y dx trapezoidal','dy/dx analytical','\int y dx analytical')
legend('location','northwest')
ylabel('see legend')
xlabel('x')

function [dydx,xMP,F]=ApproxCalculus(x,y,x0,F0)
%Differentiate using finite difference
dydx=diff(y,1,2)./diff(x,1,2);
%Calculate location of points where dydx applies
xMP=(x(1,1:end-1)+x(1,2:end))/2;
%Integrate using trapezoidal rule
F=cumtrapz(x,y,2);
%Apply constraint
F=F-F(x==x0)+F0;