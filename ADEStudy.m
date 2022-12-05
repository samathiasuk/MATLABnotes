function ADEStudy
%A sequence of different subfunctions to explore stability and
%numerical diffusion associated with finite difference solutions of
%the advection diffusion equation.

%Define model parameters
L=2; %m
a=1; %m/day
b=0.01; %m^2/day
uI=0; %mg/l
u0=1; %mg/l

%Define location of points for plotting
x=linspace(0,L,50);
%Define values of time to be used for plotting
t=[1];

%Evaluate the exact solution
uA=AnalSol(x,t,a,b,uI,u0);

%Specify a Peclet number and Courant number
Pe=2;
Cr=1/(2+Pe);
Cr=0.5;

%Evaluate the explicit time-stepping with central differences in space scheme
uFD(:,1)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ETCS');
%Evaluate the implicit time-stepping with central differences in space scheme
uFD(:,2)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ITCS');
%Evaluate the implicit time-stepping with backward differences in space scheme
uFD(:,3)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ITBS');
%Evaluate the explicit time-stepping with backward differences in space scheme
%uFD(:,4)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ETBS');
%Evaluate the ODE45 with central differences in space scheme
uFD(:,5)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ODE45CS');
%Evaluate the ODE45 with backward differences in space scheme
uFD(:,6)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ODE45BS');

%Evaluate the associated analytical solutions with numerical diffusion
uAcorr(:,1)=AnalSol(x,t,a,b*(1-Pe^2*Cr/2),uI,u0);
uAcorr(:,2)=AnalSol(x,t,a,b*(1+Pe^2*Cr/2),uI,u0);
uAcorr(:,3)=AnalSol(x,t,a,b*(1+Pe*(1+Pe*Cr)/2),uI,u0);
uAcorr(:,4)=AnalSol(x,t,a,b*(1+Pe*(1-Pe*Cr)/2),uI,u0);
uAcorr(:,5)=AnalSol(x,t,a,b,uI,u0);
uAcorr(:,6)=AnalSol(x,t,a,b*(1+Pe/2),uI,u0);

%Plot results
figure(1)
hold off
plot(x,uA,'ko')
hold on
plot(x,uFD)
plot(x,uAcorr,'.')
xlabel('Distance (m)')
ylabel('Concentration (mg/l)')
title('Concentration profile after one day')
legend('Analytical solution','ETCS','ITCS','ITBS','ETBS','ODE45 with CS','ODE45 with BS')

%**********************************************************************

function u=AnalSol(x,t,a,b,uI,u0)
%Subfunction containing the analytical solution
[x,t]=ndgrid(x,t);
Pe=a*x/b;
F1=erfc((x-a*t)./2./sqrt(b*t))/2;
F2=exp(Pe).*erfc((x+a*t)./2./sqrt(b*t))/2;
F=F1+F2;
ind=Pe>700;
F(ind)=F1(ind);
u=(u0-uI)*F+uI;

%**********************************************************************

function u=FDSol(xI,tI,a,b,uI,u0,Pe,Cr,Scheme)
%Subfunction containing the finite difference solutions

%Determine the space-step and time-step using a specified
%Courant and Peclet number
dx=Pe*b/a;
dt=Cr*dx^2/b;

%Determine number of nodes to solve for
Nx=ceil((xI(end)-xI(1))/dx)+1;
Nt=ceil(tI(end)/dt)+1;

%Determine the location of solution points in x and t
x=[0:dx:dx*(Nx-1)]+xI(1);
t=[0:dt:dt*(Nt-1)];

%Initialise the solution matrix
u=zeros(Nx,Nt);

%Apply initial condition
u(:,1)=uI;

%Apply boundary conditions
u(1,:)=u0;
u(Nx,:)=uI;

switch Scheme
    case 'ODE45CS'
        %Use ODE45 with central differences
        [t,u]=ode45(@odefunODE45CS,[0 tI],u(:,1),[],a,b,dx);
        %Note that u will need to be transposed to be the right way around
        u=u';
    case 'ODE45BS'
        %Use ODE45 with backward differences
        [t,u]=ode45(@odefunODE45BS,[0 tI],u(:,1),[],a,b,dx);
        %Note that u will need to be transposed to be the right way around
        u=u';
    otherwise
        %Determine A, B and C coefficients
        switch Scheme
            case 'ETCS'
                %Explicit time-stepping central difference in space
                A=a*dt/2/dx+b*dt/dx^2;
                B=1-2*b*dt/dx^2;
                C=-a*dt/2/dx+b*dt/dx^2;
            case 'ITCS'
                %Implicit time-stepping central difference in space
                A=-a*dt/2/dx-b*dt/dx^2;
                B=1+2*b*dt/dx^2;
                C=a*dt/2/dx-b*dt/dx^2;
            case 'ITBS'
                %Implicit time-stepping backward difference in space
                A=-a*dt/dx-b*dt/dx^2;
                B=1+a*dt/dx+2*b*dt/dx^2;
                C=-b*dt/dx^2;
            case 'ETBS'
                %Explicit time-stepping backward difference in space
                A=a*dt/dx+b*dt/dx^2;
                B=1-a*dt/dx-2*b*dt/dx^2;
                C=b*dt/dx^2;
        end
        
        %Determine Jacobian matrix
        J=spdiags(ones(Nx,1)*[A B C],[-1 0 1],Nx,Nx);
        
        %Apply boundary conditions
        J(1,:)=0;
        J(1,1)=1;
        J(Nx,:)=0;
        J(Nx,Nx)=1;
        
        %Solve problem
        switch Scheme(1:2)
            case 'ET'
                %Explicit time-stepping
                for n=1:Nt-1
                    u(:,n+1)=J*u(:,n);
                end
            case 'IT'
                %Implicit time-stepping
                for n=1:Nt-1
                    u(:,n+1)=J\u(:,n);
                end
        end
end

%Interpolate solution to plotting points
u=interp2(t,x',u,tI,xI');

%**********************************************************************

function dudt=odefunODE45CS(t,u,a,b,dx)
%ODE function for the central differencing scheme
dudt=zeros(size(u));
i=2:numel(u)-1;
dudt(i)=-a*(u(i+1)-u(i-1))/2/dx+b*(u(i-1)-2*u(i)+u(i+1))/dx^2;

%**********************************************************************

function dudt=odefunODE45BS(t,u,a,b,dx)
%ODE function for the backward differencing scheme
dudt=zeros(size(u));
i=2:numel(u)-1;
dudt(i)=-a*(u(i)-u(i-1))/dx+b*(u(i-1)-2*u(i)+u(i+1))/dx^2;
