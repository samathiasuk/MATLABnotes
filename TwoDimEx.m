function TwoDimEx
%Solves the two-dimensional diffusion of iron-oxide from clasts in a
%conglomerate using ODE15s
%
%Lx (cm) - Length of domain in the x-direction
%Ly (cm) - Length of domain in the x-direction
%rMU (cm) - Mean clast radius
%rSD (cm) - Standard deviation of clast radius
%NC (-) - Number of clasts studied
%DaM (cm^2/s) - Apparent diffusion coefficient in the rock matrix
%DaC (cm^2/s) - Apparent diffusion coefficient in the clasts
%Da (cm^2/s) - Apparent diffusion coefficient
%Dx (cm^2/s) - Apparent diffusion coefficient in the x-direction
%Dy (cm^2/s) - Apparent diffusion coefficient in the y-direction
%C0M (-) - Initial iron oxide mass-fraction in the rock matrix
%C0C (-) - Initial iron oxide mass-fraction in the clasts
%C0 (-) - Initial iron oxide mass-fraction
%Nx (-) - Number of points to solve in the x-direction
%x (m) - location of points being solved for in x-direction
%xB (m) - location of finite difference block boundaries in x-direction
%Ny (-) - Number of points to solve in the y-direction
%y (m) - location of points being solved for in y-direction
%yB (m) - location of finite difference block boundaries in y-direction
%ClastMap (-) - A logic array stating whether or not a node is within a clast
%tkYr (1000 years) - Time in 1000s of years
%t (s) - Time in seconds

%Define model parameters
Lx=100;
Ly=100;
rMU=3;
rSD=5;
NC=40;
DaM=1e-7;
DaC=1e-9;
C0M=0.01;
C0C=0.10;

%Define spatial grid
Nx=150;
Ny=155;
xB=linspace(0,Lx,Nx+1);
yB=linspace(0,Ly,Ny+1);
x=(xB(2:Nx+1)+xB(1:Nx))/2;
y=(yB(2:Ny+1)+yB(1:Ny))/2;

[xB,~]=ndgrid(xB,y);
[~,yB]=ndgrid(x,yB);
[x,y]=ndgrid(x,y);

%Generate clasts

%Initialise the ClastMap array
ClastMap=zeros(size(x));
%Start the clast counter
ClastNo=0;
while ClastNo<NC
    %Determine centre of clast
    xC=rand(1,1)*Lx;
    yC=rand(1,1)*Ly;
    %Determine radius of clast
    rC=rMU+rSD*abs(randn(1,1));
    %Calculate radial distance from point
    r=sqrt((xC-x).^2+(yC-y).^2);
    %Update the clast map check array
    ClastMapCheck=ClastMap;
    %Mark clast in clast map
    ClastMapCheck(r<rC)=ClastMapCheck(r<rC)+1;
    if max(ClastMapCheck(:))<=1
        %New clast does not overlap existing clasts
        ClastMap=ClastMapCheck;
        %Update clast counter
        ClastNo=ClastNo+1;
    end
end

%Define initial condition
C0=zeros(size(x))+C0M;
C0(ClastMap==1)=C0C;
%Define diffusion coeffients
Da=zeros(size(x))+DaM;
Da(ClastMap==1)=DaC;

%Interpolate the diffusive fluxes to the xB values
Dax=interp1(x(:,1),Da,xB(:,1));
%Interpolate the diffusive fluxes to the yB values
Day=interp1(y(1,:),Da',yB(1,:))';

%Define times to solve for in years
tYR=[0 1 10 100];
%Convert times to seconds
t=tYR*365*24*60^2;

%Vectorize the initial condition
C0=reshape(C0,Nx*Ny,1);

% %Ignore options for the ODE solver
% options=[];

%Define and set the Jacobian pattern
JPat=spdiags(ones(Nx*Ny,5),[-Nx -1 0 1 Nx],Nx*Ny,Nx*Ny);
options=odeset('JPattern',JPat);

%Solve problem
h=waitbar(0,'Please wait...');
[t,C]=ode15s(@MYodefun,t,C0,options,x,xB,y,yB,Nx,Ny,Dax,Day,h,t(end));
close(h)

%Plot the results
figure(1)
for n=1:4
    subplot(2,2,n)
    %Show results as a coloured surface
    surf(x,y,reshape(C(n,:)',Nx,Ny))
    %Set the x and y axis limits
    axis([0 Lx 0 Ly])
    %Set the color-scale limites
    caxis([C0M C0C])
    %Display the color bar
    colorbar
    %Use linear interpolation in the colour shading
    shading interp
    %Show a plan view of the surface
    view(2)
    %Label x and y axis
    xlabel('Distance (cm)')
    ylabel('Distance (cm)')
    %Add a title to each subplot
    title(['Iron oxide mass fraction after ' num2str(tYR(n)) ' years'])
end


function dCdt=MYodefun(t,C,x,xB,y,yB,Nx,Ny,Dax,Day,h,tmax)

%De-vectorize the C vector
C=reshape(C,Nx,Ny);

%Calculate the diffusive fluxes
qx=-Dax(2:Nx,:).*diff(C,1,1)./diff(x,1,1);
qy=-Day(:,2:Ny).*diff(C,1,2)./diff(y,1,2);

%Apply zero flux boundary conditions
qx=[zeros(1,Ny); qx; zeros(1,Ny)];
qy=[zeros(Nx,1) qy zeros(Nx,1)];

%Apply the mass conservation equation
dCdt=-diff(qx,1,1)./diff(xB,1,1)-diff(qy,1,2)./diff(yB,1,2);

%Vectorize the C vector
dCdt=reshape(dCdt,Nx*Ny,1);

%Monitor progress
waitbar(t/tmax,h)


