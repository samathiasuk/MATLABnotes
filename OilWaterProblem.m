function OilWaterProblem
%Solution for water injection in an oil reservoir problem


%Define model parameters
aw=2;
an=2;
b=0.5;
pc0=20000;
muw=1e-3;
mun=2e-3;
k=1e-13;
qt=0.2;
phi=0.2;
L=120;
SwI=0.2;

%Define spatial grid
N=150;
xB=linspace(0,L,N+1)';
x=(xB(1:N,1)+xB(2:N+1,1))/2;

%Define the times of interest
t=[0 20 40 60 80];
%Convert t to seconds
t=t*24*60^2;

%Convert qt to seconds
qt=qt/24/60^2;

%Define intial condition vector for the ode solver
SwI_vec=zeros(size(x))+SwI;

%Define and set the Jacobian pattern
JPat=spdiags(ones(N,3),[-1 0 1],N,N);
options=odeset('JPattern',JPat);

%Apply ode15s to obtain the finite difference solution
[t,Sw]=ode15s(@MYodefun,t,SwI_vec,options,x,xB,aw,an,b,pc0,muw,mun,k,qt,phi);

[xA,SwA]=AnalyticalSol(0,0,an,aw,1,1,mun,muw,SwI,qt,phi,t);

%Plot results and compare with the analytical solution
figure(1)
hold off
plot(x,Sw)
hold on
plot(xA,SwA,'--')
xlim([0 L])
ylabel('Water saturation')
xlabel('Distance from injection (m)')
legend('0 days','20 days','40 days','60 days','80 days','Analytical solution')

%**************************************************************************

function dSwdt=MYodefun(t,Sw,x,xB,aw,an,b,pc0,muw,mun,k,qt,phi)
%Calculate derivatives
dSwdx=diff(Sw,1,1)./diff(x,1,1);
%Interpolate Sw to the internal xB points
Sw=interp1(x,Sw,xB(2:end-1,1));
%Calculate auxillary variables
krw=Sw.^aw;
Sn=1-Sw;
krn=Sn.^an;
pc=pc0.*(Sw.^(-1/b)-1).^(1-b);
dpcdSw=pc*(1-b)/b./Sw./(Sw.^(1/b)-1);
fw=1./(1+krn*muw./krw/mun);
%Calculate capillary flux at internal xB points
qc=-k*krw/muw.*dpcdSw.*dSwdx;
%Calculate water flux at internal xB points
qw=qt*fw-(1-fw).*qc;
%Apply boundary conditions
qw=[qt;qw;qt*fw(end)];
%Calculate the flux divergence
dqwdx=diff(qw,1,1)./diff(xB,1,1);
%Solve the mass conservation statement
dSwdt=-dqwdx/phi;

%**************************************************************************
%*******************    Analytical Solution Stuff    **********************
%**************************************************************************

function [x,Sg]=AnalyticalSol(Sar,Sgc,m,n,kra0,krg0,mua,mug,Sg0,q0,phi,t)
%Find the saturation at the shock front
Sg1=fminsearch(@ShockFun,0.9*(1-Sar),[],Sar,Sgc,m,n,kra0,krg0,mua,mug,Sg0);
%generate values from the shock front to the boundary value
Sg=linspace(Sg1,1-Sar,100);

%Generate the self-similar solution
[dfgdSg,fg]=PowerLawFun(Sg,Sar,Sgc,m,n,kra0,krg0,mua,mug);
[~,fg0]=PowerLawFun(Sg0,Sar,Sgc,m,n,kra0,krg0,mua,mug);

%Add Sg0 to the vectors
Sg=[Sg0 Sg0 Sg];
fg=[fg0 fg0 fg];
dfgdSg=[dfgdSg(1)*10 dfgdSg(1) dfgdSg];

%Calculate locations for different times
z=dfgdSg;
[t,z]=ndgrid(t,z);

%distance in m
x=z*q0.*t/phi; 

%**************************************************************************

function [F,dfgdSg,fg]=ShockFun(Sg,Sar,Sgc,m,n,kra0,krg0,mua,mug,Sg0)
[dfgdSg0,fg0]=PowerLawFun(Sg0,Sar,Sgc,m,n,kra0,krg0,mua,mug);
[dfgdSg,fg,kra,krg]=PowerLawFun(Sg,Sar,Sgc,m,n,kra0,krg0,mua,mug);
F=abs((fg-fg0)./(Sg-Sg0)./dfgdSg-1);

%**************************************************************************

function [dfgdSg,fg,kra,krg]=PowerLawFun(Sg,Sar,Sgc,m,n,kra0,krg0,mua,mug)
%See Mathias et al. (2011, WRR)

%Eqs. 42 to 43 of paper
kra=kra0*((1-Sg-Sar)/(1-Sgc-Sar)).^m;
krg=krg0*((Sg-Sgc)/(1-Sgc-Sar)).^n;
%eqn 41 of paper
fg=(1+(mug*kra)./(mua*krg)).^-1;
dfgdSg=fg.*(1-fg).*(n.*(1-Sg-Sar)+m.*(Sg-Sgc))./(Sg-Sgc)./(1-Sg-Sar);

%Apply Sgc constraint
ind=Sg<=Sgc;
kra(ind)=kra0;
krg(ind)=0;
fg(ind)=0;
dfgdSg(ind)=0;