function ADEproblem
%Solves the one-dimensional advection dispersion equation using a MATLAB
%ODE solver
%
%t (days) - time after start of spill
%x (km) - distance from spill
%v (km/day) - velocity of river flow
%DL (km^2/day) - longitudinal dispersion coefficient
%CI (mg/l) - initial concenctration
%Cdump (mg/l) - concentration at spill site during spill incident
%tdump (days) - duration of spill incident

%Define model parameters
CI=0.1;
Cdump=120;
tdump=1;
v=1.2;
L=10;
DL=0.23/20;

%Define spatial grid
N=50;
xB=linspace(0,L,N+1)';
x=(xB(1:N,1)+xB(2:N+1,1))/2;

%Define the times of interest
t=[0 2 4 6];

%Define intial condition vector for the ode solver
CI_vec=zeros(size(x))+CI;

%Define and set the Jacobian pattern
JPat=spdiags(ones(N,3),[-1 0 1],N,N);
options=odeset('JPattern',JPat);

%Apply ode15s to obtain the finite difference solution
[t,CFD]=ode15s(@MYodefun,t,CI_vec,options,x,xB,v,DL,CI,Cdump,tdump);

%Evaluate analytical solution
C=AnalyticalSol(t,x,v,DL,CI,Cdump,tdump);

%Plot results and compare with the analytical solution
figure(1)
hold off
plot(x,C)
hold on
plot(x,CFD,'.--')
ylabel('Product concentration (mg/l)')
xlabel('Distance from spill (km)')
legend('0 days (analytical)','2 days (analytical)',...
    '4 days (analytical)','6 days (analytical)'...
    ,'0 days (finite difference)','2 days (finite difference)'...
    ,'4 days (finite difference)','6 days (finite difference)')

%**************************************************************************

function dCdt=MYodefun(t,C,x,xB,v,DL,CI,Cdump,tdump)
%Determine C0
C0=Cdump;
C0(t>tdump)=CI;
%Calculate derivatives and include C=C0 at x=0 boundary condition
dCdx=diff([C0;C],1,1)./diff([xB(1);x],1,1);
%Interpolate C to the xB points
CB=interp1(x,C,xB);
%Impose fixed concentration boundary
CB(1)=C0;
%Impose zero gradient boundary
CB(end)=C(end);
%Calculate Dispersive fluxes and include q=0 at x=L boundary condition
%Also add on the advective flux, v*CB
q=[-DL*dCdx;0]+v*CB;
%Calculate the flux divergence
dqdx=diff(q,1,1)./diff(xB,1,1);
%Solve the mass conservation statement
dCdt=-dqdx;

%**************************************************************************

function C=AnalyticalSol(t,x,v,DL,CI,Cdump,tdump)
%Analytical solution for the advection dispersion problem
%t (days) - time after start of spill
%x (km) - distance from spill
%v (km/day) - velocity of river flow
%DL (km^2/day) - longitudinal dispersion coefficient
%CI (mg/l) - initial concenctration
%Cdump (mg/l) - concentration at spill site during spill incident
%tdump (days) - duration of spill incident
[t,x]=ndgrid(t,x);
F1=FFun(t*v./x,v*x/DL);
F2=FFun((t-tdump)*v./x,v*x/DL);
F2(t<tdump)=0;
F=F1-F2;
C=(Cdump-CI)*F+CI;

%**************************************************************************

function F=FFun(z,P)
z(z<0)=0;
F=(erfc(sqrt(P/4./z).*(1-z))+exp(P).*erfc(sqrt(P/4./z).*(1+z)))/2;
ind=P>700;
F(ind)=erfc(sqrt(P(ind)/4./z(ind)).*(1-z(ind)))/2;