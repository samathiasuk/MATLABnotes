function NNESMOwaves



%Define number of solution points
nx=500;
ny=100;
nyi=40;
nxs=160;
nys=2;
nt=1800;
dt=0.001;
dx=10;
dy=10;
x=[dx:dx:nx*dx]';
y=[dy:dy:ny*dy];
t=dt:dt:nt*dt;

%Define model parameters
muA=1e10;
rhoA=2500;
muB=2.5e10;
rhoB=2777;
mu=muA+zeros(nx,ny);
rho=rhoA+zeros(nx,ny);
mu(:,1:nyi)=muB;
rho(:,1:nyi)=rhoB;

mux=1./(1./mu(1:end-1,:)+1./mu(2:end,:));
muy=1./(1./mu(:,1:end-1)+1./mu(:,2:end));
mux=[mu(1,:); mux];
muy=[mu(:,1) muy];
rhox=1./(1./rho(1:end-1,:)+1./rho(2:end,:));
rhoy=1./(1./rho(:,1:end-1)+1./rho(:,2:end));
rhox=[rho(1,:); rhox];
rhoy=[rho(:,1) rhoy];





%Define source term
tau=0.02;
t0=2*tau;
v0=exp(-((t-t0)/tau).^2)/tau;




%Initialise v and s arrays
v=zeros(nx,ny,2);
v(nxs,nys)=v0(1);
sx=zeros(nx,ny,2);
sy=zeros(nx,ny,2);


%Solve problem
%Set n=1 so that the solution always updates itself
n=1;
for m=1:nt-1
    %From slide 21
    i=2:nx;
    j=2:ny;
    dsxdx=(sx(i,j,n)-sx(i-1,j,n))/dx;
    dsydy=(sy(i,j,n)-sy(i,j-1,n))/dy;
    dvdt=(dsxdx+dsydy)./rho(i,j);
    
    v(i,j,n+1)=v(i,j,n)+dvdt*dt;
    %Apply source at nxs and nys
    v(nxs,nys,n+1)=v(nxs,nys,n+1)+v0(m+1);
    %Apply v=0 at x=0 and y=0
    v(1,:,n+1)=0;
    v(:,1,n+1)=0;
    
    %From slide 21
    i=1:nx-1;
    j=1:ny-1;
    dvdx=(v(i+1,j,n+1)-v(i,j,n+1))/dx;
    dvdy=(v(i,j+1,n+1)-v(i,j,n+1))/dy;
    dsxdt=mux(i,j).*dvdx;
    dsydt=muy(i,j).*dvdy;
    
    sx(i,j,n+1)=sx(i,j,n)+dsxdt*dt;
    sy(i,j,n+1)=sy(i,j,n)+dsydt*dt;
    
    %Apply s=0 at x=Lx
    sx(nx,:,n+1)=0;
    %Apply s=0 at y=Ly
    sy(:,ny,n+1)=0;
    
    
    %Update solution
    sx(:,:,n)=sx(:,:,n+1);
    sy(:,:,n)=sy(:,:,n+1);
    v(:,:,n)=v(:,:,n+1);

    mPLOT=3;
    if round(m/mPLOT)*mPLOT==m
        %Now plot results
        figure(1)
        clf
        surf(x,y,v(:,:,n)','EdgeColor','none')
        
        view(2)
        shading interp
        caxis([-1 1])
        title(t(m))
        drawnow
    end
end
