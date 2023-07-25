function MATLABsession6_Assignment
%This file contains instructions requested by MATLAB Notes - Session 6

%Compute instructions associated with the groundwater flow example
GroundwaterExample

%Compute instructions associated with the chemical transport example
ChemicalTransportExample

%**************************************************************************

function GroundwaterExample
%Solves the one-dimensional confined aquifer flow equation using a MATLAB
%ODE solver
%T (m^2/day) - Transmissivity
%S (-) - Storativity
%L (m) - Length of the aquifer
%hI (m) - Initial hydraulic head
%h0 (m) - Boundary head
%x (m) - location of points being solved for
%xB (m) - location of finite difference block boundaries
%t (days) - time
%Define model parameters
T=800;
S=0.01;
L=1600;
hI=50;
h0=53;
%Define spatial grid
N=20;
xB=[0;logspace(-3,0,N)'*L];
x=(xB(1:N,1)+xB(2:N+1,1))/2;
%Define the times of interest
t=[0 0.1 1 10 100];
%Define intial condition vector for the ode solver
hI_vec=zeros(size(x))+hI;
%Define and set the Jacobian pattern
JPat=spdiags(ones(N,3),[-1 0 1],N,N);
options=odeset('JPattern',JPat);
%Apply ode15s to obtain the finite difference solution
[t,hFD]=ode15s(@GroundwaterODEfun,t,hI_vec,options,x,xB,S,T,h0);
%Evaluate analytical solution
h=GroundwaterAnalyticalSol(t,x,S,T,hI,h0);
%Plot results and compare with the analytical solution
figure(1)
clf
hold on
plot(x,h)
%Reset colour order for plotting
set(gca,'ColorOrderIndex',1)
plot(x,hFD,'o--')
xlabel('Distance (m)')
ylabel('Hydraulic head (m)')
legend('0 days (analytical)','0.1 days (analytical)',...
    '1 days (analytical)','10 days (analytical)','100 days (analytical)'...
    ,'0 days (finite difference)','0.1 days (finite difference)'...
    ,'1 days (finite difference)','10 days (finite difference)'...
    ,'100 days (finite difference)')
%Make a nice box around the graphs
box on

%**************************************************************************

function dhdt=GroundwaterODEfun(t,h,x,xB,S,T,h0)
%Calculate derivatives and include h=h0 at x=0 boundary condition
dhdx=diff([h0;h],1,1)./diff([xB(1);x],1,1);
%Calculate Darcy fluxes and include Q=0 at x=L boundary condition
Q=[-T*dhdx;0];
%Calculate the flux divergence
dQdx=diff(Q,1,1)./diff(xB,1,1);
%Solve the mass conservation statement
dhdt=-dQdx/S;

%**************************************************************************

function h=GroundwaterAnalyticalSol(t,x,S,T,hI,h0)
[t,x]=ndgrid(t,x);
h=(h0-hI)*erfc(sqrt(S*x.^2/4/T./t))+hI;

%**************************************************************************
%**************************************************************************
%**************************************************************************

function ChemicalTransportExample
%Solves the one-dimensional advection dispersion equation using a MATLAB
%ODE solver
%
%t (days) - time after start of spill
%x (km) - distance from spill
%v (km/day) - velocity of river flow
%DL (km^2/day) - longitudinal dispersion coefficient
%CI (mg/l) - initial concenctration
%C00 (mg/l) - concentration at spill site during spill incident
%t0 (days) - duration of spill incident

%Define model parameters
CI=0.1;
C00=120;
t0=1;
v=1.2;
L=10;
DL=0.23;

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
[t,CFD]=ode15s(@ChemicalTransportODEfun,t,CI_vec,options,x,xB,v,DL,CI,C00,t0);

%Evaluate analytical solution
C=ChemicalTransportAnalyticalSol(t,x,v,DL,CI,C00,t0);

%Plot results and compare with the analytical solution
figure(2)
clf
hold on
plot(x,C)
%Reset colour order for plotting
set(gca,'ColorOrderIndex',1)
plot(x,CFD,'o--')
%Make a nice box around the graphs
box on
ylabel('Product concentration (mg l^-^1)')
xlabel('Distance from spill (km)')
legend('0 days (analytical)','2 days (analytical)',...
    '4 days (analytical)','6 days (analytical)'...
    ,'0 days (finite difference)','2 days (finite difference)'...
    ,'4 days (finite difference)','6 days (finite difference)')

%**************************************************************************

function dCdt=ChemicalTransportODEfun(t,C,x,xB,v,DL,CI,C00,t0)
%Determine C0
C0=C00;
C0(t>t0)=CI;
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

function C=ChemicalTransportAnalyticalSol(t,x,v,DL,CI,C00,t0)
%Analytical solution for the advection dispersion problem
%t (days) - time after start of spill
%x (km) - distance from spill
%v (km/day) - velocity of river flow
%DL (km^2/day) - longitudinal dispersion coefficient
%CI (mg/l) - initial concenctration
%C00 (mg/l) - concentration at spill site during spill incident
%t0 (days) - duration of spill incident
[t,x]=ndgrid(t,x);
F1=FFun(t*v./x,v*x/DL);
F2=FFun((t-t0)*v./x,v*x/DL);
F2(t<t0)=0;
F=F1-F2;
C=(C00-CI)*F+CI;

%**************************************************************************
function F=FFun(z,P)
z(z<0)=0;
F=(erfc(sqrt(P/4./z).*(1-z))+exp(P).*erfc(sqrt(P/4./z).*(1+z)))/2;
ind=P>700;
F(ind)=erfc(sqrt(P(ind)/4./z(ind)).*(1-z(ind)))/2;

