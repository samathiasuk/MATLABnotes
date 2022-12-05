function BurgersStudy
%A sequence of different subfunctions to solve the Burgers' equation.

%Define model parameters
nu=0.1; %m^2/day
u0=1.2; %m/day
x0=0.7; %m

%Define location of points for plotting
x=linspace(-4,6,100);
%Define values of time (days) to be used for plotting
t=[0.5 1 2 4 8];

%Evaluate the exact solution
uA=AnalSol(x,t,nu,u0,x0);

%Evaluate the finite difference solutions
dx=0.02;
dt=0.25;

%MAXITER is the a variable containing the maximum number of iterations

tic
%Apply semi-implicit method by setting MAXITER to 1
uSITBS=FDSol(x,t,nu,x0,u0,dx,dt,1);
toc

tic
%Apply implicit method by setting MAXITER to 100
uITBS=FDSol(x,t,nu,x0,u0,dx,dt,100);
toc

tic
%Apply ODE15s with backward differencing by setting MAXITER to []
uODE15sBACK=FDSol(x,t,nu,x0,u0,dx,dt,[],'backwards');
toc

tic
%Apply ODE15s with central differencing by setting MAXITER to NaN
uODE15sCENT=FDSol(x,t,nu,x0,u0,dx,dt,[],'central');
toc

%Plot results
figure(1)
hold off
plot(x,uA,'.')
hold on
%You can use line types in addition to color in your legend
%if you plot the first time of each model type on its own first
plot(x,uSITBS(:,1),'-')
plot(x,uITBS(:,1),'--')
plot(x,uODE15sBACK(:,1),':')
plot(x,uODE15sCENT(:,1),'-.')
%Now plot everything else, but these lines won't have a legend entry
plot(x,uSITBS,'-')
plot(x,uITBS,'--')
plot(x,uODE15sBACK,':')
plot(x,uODE15sCENT,'-.')
xlabel('Distance (m)')
ylabel('Velocity (m/day)')
for n=1:numel(t)
    MYlegend{n}=[num2str(t(n)) ' days'];
end
MYlegend{6}='SITBS';
MYlegend{7}='ITBS';
MYlegend{8}='ODE15s with backwards';
MYlegend{9}='ODE15s with central';
legend(MYlegend,'location','northwest')

%**********************************************************************

function u=AnalSol(x,t,nu,u0,x0)
%Subfunction containing the analytical solution
[x,t]=ndgrid(x,t);
x1=(x+x0-u0*t)./sqrt(4*nu*t);
x2=(x-x0-u0*t)./sqrt(4*nu*t);
x3=(x0+x)./sqrt(4*nu*t);
x4=(x0-x)./sqrt(4*nu*t);
c=u0/2/nu;
T3=exp(c*(x+x0-u0*t/2));
T4=exp(c*(x-x0-u0*t/2));
u=u0./(1+(T3.*erfc(x3)+T4.*erfc(x4))./(erfc(x2)-erfc(x1)));

%**********************************************************************

function u=FDSol(xI,tI,nu,x0,u0,dx,dt,MAXITER,Scheme)
%Subfunction containing the finite difference solutions

%Note that it is possible to overload the input arguments

%Determine number of nodes to solve for
Nx=ceil((xI(end)-xI(1))/dx)+1;
Nt=ceil(tI(end)/dt)+1;

%Determine the location of solution points in x and t
x=[0:dx:dx*(Nx-1)]+xI(1);
t=[0:dt:dt*(Nt-1)];

%Define initial condition
uI=ones(size(x'))*u0;
uI(abs(x)>x0)=0;

if isempty(MAXITER)
    %Define Jacobian pattern for ode15s
    options=odeset('JPattern',spdiags(ones(Nx,3),[-1 0 1],Nx,Nx));
    switch Scheme
        case 'backwards'
            %Use ODE15s with backwards differences
            [t,u]=ode15s(@odefunBACK,[0 tI],uI,options,nu,dx);
        case 'central'
            %Use ODE15s with central differences
            [t,u]=ode15s(@odefunCENT,[0 tI],uI,options,nu,dx);
    end
    %Note that u will need to be transposed to be the right way around
    u=u';
else
    %Initialise the solution matrix
    u=zeros(Nx,Nt);
    
    %Apply initial condition
    u(:,1)=uI;
    
    %Apply boundary conditions
    u(1,:)=0;
    u(Nx,:)=0;
    
    %Define convergence tolerance
    MyTol=1e-6;
    
    %Step through time-steps
    for n=1:Nt-1
        
        %Set first guess
        uold=u(:,n);
        
        %Reset iteration counter
        m=1;
        Err=inf;
        
        %Step through Newton iterations
        while m<=MAXITER&Err>MyTol
            %Implicit time-stepping backward difference in space
            A=-uold*dt/dx-nu*dt/dx^2;
            B=1+uold*dt/dx+2*nu*dt/dx^2;
            C=-ones(Nx,1)*nu*dt/dx^2;
            
            %Now that A and C are variables we need to be more careful
            %about how spdiags applies these.
            %Consequently we need to move them around like this:
            A=[A(2:end);0];
            C=[0;C(1:end-1)];
            %Determine Jacobian matrix
            J=spdiags([A B C],[-1 0 1],Nx,Nx);
            
            %Apply boundary conditions
            J(1,:)=0;
            J(1,1)=1;
            J(Nx,:)=0;
            J(Nx,Nx)=1;
            
            %Solve problem
            u(:,n+1)=J\u(:,n);
            
            %Check for convergence
            Err=abs(1-uold./u(:,n+1));
            %Exlcude Nan and inf
            Err(isnan(Err)|isinf(Err))=[];
            %Select maximum error
            Err=max(Err);
            
            %Update guess
            uold=u(:,n+1);
            %Update iteration count
            m=m+1;
            
            %Report if not converged
            if m==MAXITER
                disp(['Not converged at t = ' num2str(t(n+1)) ' days'])
            end
            
        end
    end
end

%Interpolate solution to plotting points
u=interp2(t,x',u,tI,xI');


%**********************************************************************

function dudt=odefunBACK(t,u,nu,dx)
%ODE function for the backward differencing scheme
dudt=zeros(size(u));
i=2:numel(u)-1;
dudt(i)=-u(i).*(u(i)-u(i-1))/dx+nu*(u(i-1)-2*u(i)+u(i+1))/dx^2;

%**********************************************************************

function dudt=odefunCENT(t,u,nu,dx)
%ODE function for the central differencing scheme
dudt=zeros(size(u));
i=2:numel(u)-1;
dudt(i)=-u(i).*(u(i+1)-u(i-1))/2/dx+nu*(u(i-1)-2*u(i)+u(i+1))/dx^2;
