function OneDimEx
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
[t,hFD]=ode15s(@MYodefun,t,hI_vec,options,x,xB,S,T,h0);

%Evaluate analytical solution
h=AnalyticalSol(t,x,S,T,hI,h0);




%Plot results and compare with the analytical solution
figure(1)
hold off
plot(x,h)
hold on
plot(x,hFD,'.--')
xlabel('Distance (m)')
ylabel('Hydraulic head (m)')
legend('0 days (analytical)','0.1 days (analytical)',...
    '1 days (analytical)','10 days (analytical)','100 days (analytical)'...
    ,'0 days (finite difference)','0.1 days (finite difference)'...
    ,'1 days (finite difference)','10 days (finite difference)'...
    ,'100 days (finite difference)')

%**************************************************************************

function dhdt=MYodefun(t,h,x,xB,S,T,h0)
%Calculate derivatives and include h=h0 at x=0 boundary condition
dhdx=diff([h0;h],1,1)./diff([xB(1);x],1,1);
%Calculate Darcy fluxes and include Q=0 at x=L boundary condition
Q=[-T*dhdx;0];
%Calculate the flux divergence
dQdx=diff(Q,1,1)./diff(xB,1,1);
%Solve the mass conservation statement
dhdt=-dQdx/S;

%**************************************************************************

function h=AnalyticalSol(t,x,S,T,hI,h0)
[t,x]=ndgrid(t,x);
h=(h0-hI)*erfc(sqrt(S*x.^2/4/T./t))+hI;