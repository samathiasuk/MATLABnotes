%% NNESMO: Session 7 - Solving two dimensional PDEs using ODE solvers
%
% Dr Simon A Mathias
%
% Department of Earth Sciences
%
% Durham University
%
%% This session's learning objectives
%
% At the end of this session you should be able to:
%
% * Use |ndgrid| to define the spatial coordinates for n-dimensional
% finite difference problems.
% * Use |reshape| to vectorize n-dimensional problems.
% * Use |ode45| and |ode15s| to solve n-dimensional partial differential
% equations.
% * Specify the Jacobian pattern for a two-dimensional finite difference
% problem.
% * Solve a two-dimensional diffusion problem with spatially varying
% diffusion coefficient and initial condition.
%
%% Iron oxide diffusion in a newly formed conglomerate
%
% This practial is based on an exercise previously written by Lara Kalnins
% at Durham University in 2014.
%
% In this study we will look at the diffusion of iron oxide in a
% poorly-sorted conglomerate where the larger clasts are of significantly 
% different composition as compared to the host matrix. Such a problem can be
% meaningfully investigated using a finite difference solution of the
% two-dimensional diffusion equation.
%
% The governing equation for two-dimensional diffusion takes the form
%
% $\displaystyle \frac{\partial c}{\partial t}=-\frac{\partial
% q_x}{\partial x}-\frac{\partial q_y}{\partial y}$
%
% where $c$ [-] is the mass fraction of iron oxide in the material being
% studied, $t$ [T] is time, $x$ [L] and $y$ [L] are distances in the horizontal
% and vertical direction, respectively, and the diffusive fluxes, $q_x$ $[\mathrm{LT}^{-1}]$ and
% $q_y$ $[\mathrm{LT}^{-1}]$ are found from:
%
% $\displaystyle q_x=-D_A(x,y)\frac{\partial c}{\partial x}$
%
% $\displaystyle q_y=-D_A(x,y)\frac{\partial c}{\partial y}$
%
% where $D_A$ is the apparent diffusion coefficient.
%
% The initial and boundary conditions can be written as follows:
%
% $\begin{array}{llll}
% c=c_0(x,y), & 0\leq x \leq L_x, & 0\leq y \leq L_y, & t=0\\
% q_x=0, & x=0, & 0\leq y \leq L_y, & t>0\\
% q_x=0, & x=L_x, & 0\leq y \leq L_y, & t>0\\
% q_y=0, & 0\leq x \leq L_x, & y=0, & t>0\\
% q_y=0, & 0\leq x \leq L_x, & y=L_y, & t>0\\
% \end{array}$
%
% where $c_0$ [-] is the initial distribution of iron oxide mass fraction and  
% $L_x$ [L] and $L_y$ [L] are the horizontal and vertical extents of the
% domain, which are both assumed to be 100 cm in this case. 
%
% The apparent diffusion coefficients in the matrix and clasts  are  assumed
% to be $10^{-7}$ $\mathrm{cm}^2/\mathrm{s}$ and $10^{-9}$ $\mathrm{cm}^2/\mathrm{s}$, respectively.
% The initial iron oxide mass fractions in the matrix and clasts are  assumed
% to be $0.01$ and $0.1$, respectively.
%
% There are assumed to be 40 randomly located clasts of random radii within the designated
% domain. The locations of the clast centres can be treated as  completely
% random. The radii of the clasts are assumed to be normally
% distributed with a mean of 3 cm and a standard deviation of 5 cm.
% Furthermore, all the clasts are assumed to have a radii greater than 3
% cm. It is also the case that none of the clasts overlap each other.
%
% In the previous study we used |ode45| to solve a one-dimensional partial
% differential equation (PDE). In this study will exploit a new command called 
% |reshape|, which will allow us to solve two-dimensional (or even n-dimensional) PDEs
% using |ode45|.
%
%% Spatial discretisation using finite differences
%
% Again, we are going to discretise the PDE using finite differences.
%
% Let us consider a set of $N_x$ and $N_y$ discrete points on the $x$ and
% $y$ axes, respectively:
%
% $x_i, \quad i=1,2,\dots,N_x$
%
% $y_j, \quad j=1,2,\dots,N_y$
%
% The corresponding set of mass fractions can be written as
%
% $c_{i,j}, \quad i=1,2,\dots,N_x, \quad j=1,2,\dots,N_y$
%
% The diffusive fluxes can then be written as:
%
% $\displaystyle q_{x,i+1/2,j}=-D_{A,i+1/2,j}\left(\frac{c_{i+1,j}-c_{i,j}}{x_{i+1}-x_i}\right)$
%
% $\displaystyle q_{y,i,j+1/2}=-D_{A,i,j+1/2}\left(\frac{c_{i,j+1}-c_{i,j}}{y_{j+1}-y_j}\right)$
%
% The resulting set of ODEs to be solved takes the form
%
% $\displaystyle
% \left.\frac{dc}{dt}\right|_{i,j}=-\left(\frac{q_{x,i+1/2,j}-q_{x,i-1/2,j}}{x_{i+1/2}-x_{i-1/2}}\right)
% -\left(\frac{q_{y,i,j+1/2}-q_{y,i,j-1/2}}{y_{j+1/2}-y_{j-1/2}}\right)$
%
% Also note that
%
% $x_{i+1/2}=(x_{i+1}+x_i)/2$
%
% and
%
% $y_{j+1/2}=(y_{j+1}+y_j)/2$
%
%% Development of the solution code
%
% Create a new script-file and type the following
%
% 
%   function TwoDimEx
%   %Solves the two-dimensional diffusion of iron-oxide from clasts in a
%   %conglomerate using ODE15s
%   %
%   %Lx (cm) - Length of domain in the x-direction
%   %Ly (cm) - Length of domain in the x-direction
%   %rMU (cm) - Mean clast radius
%   %rSD (cm) - Standard deviation of clast radius
%   %NC (-) - Number of clasts studied
%   %DaM (cm^2/s) - Apparent diffusion coefficient in the rock matrix
%   %DaC (cm^2/s) - Apparent diffusion coefficient in the clasts
%   %Da (cm^2/s) - Apparent diffusion coefficient
%   %Dx (cm^2/s) - Apparent diffusion coefficient in the x-direction
%   %Dy (cm^2/s) - Apparent diffusion coefficient in the y-direction
%   %C0M (-) - Initial iron oxide mass-fraction in the rock matrix
%   %C0C (-) - Initial iron oxide mass-fraction in the clasts
%   %C0 (-) - Initial iron oxide mass-fraction
%   %Nx (-) - Number of points to solve in the x-direction
%   %x (m) - location of points being solved for in x-direction
%   %xB (m) - location of finite difference block boundaries in x-direction
%   %Ny (-) - Number of points to solve in the y-direction
%   %y (m) - location of points being solved for in y-direction
%   %yB (m) - location of finite difference block boundaries in y-direction
%   %ClastMap (-) - A logic array stating whether or not a node is within a clast
%   %tkYr (1000 years) - Time in 1000s of years
%   %t (s) - Time in seconds
%   
%   %Define model parameters
%   Lx=100;
%   Ly=100;
%   rMU=3;
%   rSD=5;
%   NC=40;
%   DaM=1e-7;
%   DaC=1e-9;
%   C0M=0.01;
%   C0C=0.10;
%
% So far all we have done is written some comments and defined some
% parameters. The next step is to determine the locations on the x-axis and
% y-axis where we are going to solve for.
%
% In the above variable list we have |x|, |xB|, |y| and |yB|. 
% The |x| 
% and |y| vectors will contain $x_i$ and $y_j$, respectively. The
%  |xB| and |yB| vectors will contain the associated values for
% $x_{i+1/2}$ and $y_{j+1/2}$, respectively.
%
% Note that |xB(1)=0|, |xB(Nx+1)=Lx|, |yB(1)=0| and |yB(Ny+1)=Ly|.
%
% Add the following code to provide 100 and 105 equally spaced solution
% points in the x and y directions, respectively:
%
%   %Define spatial grid
%   Nx=100;
%   Ny=105;
%   xB=linspace(0,Lx,Nx+1);
%   yB=linspace(0,Ly,Ny+1);
%   x=(xB(2:Nx+1)+xB(1:Nx))/2;
%   y=(yB(2:Ny+1)+yB(1:Ny))/2;
%
% It is a good idea to specify slightly different numbers of points in each
% direction to avoid run-time errors in the code. If $N_x=N_y$, the code
% may still run if you have made some errors concerning direction. If
% $N_x\ne N_y$, the code should alway crash when directional errors are
% present.
%  
% It is necessary to transform these vectors to two-dimensional (2D) arrays such
% that there is a value of |xB| for every |y| value and a value of |yB| for
% every |x| value. It is then necessary to ensure there is an |x| value for
% every |y| value and vice versa. This is achieved by adding the following
% code:
%
%   [xB,~]=ndgrid(xB,y);
%   [~,yB]=ndgrid(x,yB);
%   [x,y]=ndgrid(x,y);
%
% Read the help files for |ndgrid| to find out more about how this works.
% The tilda character, |~|, can be used in function outputs to avoid
% receiving the output being replaced. This is useful here because we
% do not want to overide the vectors |x| and |y| until the third line of
% the above code.
%
% The following code can be added to develop a 2D array containing information regarding
% the locations and extents of the clasts.
%
%   %Generate clasts
%   
%   %Initialise the ClastMap array
%   ClastMap=zeros(size(x));
%   %Start the clast counter
%   ClastNo=0;
%   while ClastNo<NC
%       %Determine centre of clast
%       xC=rand(1,1)*Lx;
%       yC=rand(1,1)*Ly;
%       %Determine radius of clast
%       rC=rMU+rSD*abs(randn(1,1));
%       %Calculate radial distance from point
%       r=sqrt((xC-x).^2+(yC-y).^2);
%       %Update the clast map check array
%       ClastMapCheck=ClastMap;
%       %Mark clast in clast map
%       ClastMapCheck(r<rC)=ClastMapCheck(r<rC)+1;
%       if max(ClastMapCheck(:))<=1
%           %New clast does not overlap existing clasts
%           ClastMap=ClastMapCheck;
%           %Update clast counter
%           ClastNo=ClastNo+1;
%       end
%   end
%
% The above code uses a |while| loop to generate each clast. The |while|
% loop repeats everything contained within it until the value of |ClastNo|
% exceeds the number of required clasts, |NC|. Study the help file for
% |while| if necessary.
% 
% The first step involves randomly sampling the locations of a clast center,
% (|xC|, |yC|), using uniform random sampling from 0 to $L_x$. Read the help file for |rand|
% to find out more.
%
% Next we determine the radius of the clast, |rC|, by randomly
% sampling from the normal distribution. Read about |randn| to find out
% more. Note that here we are using only the absolute value of
% |randn(1,1)|. This ensures we only sample clasts of radii greater than
% |rMU|, as desired in the original conceptual model for this problem.
%
% In the next piece of code (i.e., |r=sqrt((xC-x).^2+(yC-y).^2|), we are 
% determining the radial distance of each
% point, $(x_i,y_j)$, from the centre of our clast using pythagoras
% theorem. Note that if |r| for a given point is |<rC|, this implies that
% the point lies within the newly generated clast.
%
% The |ClastMap| array originally contains a zero for every point,
% $(x_i,y_j)$. Zeros are another way of saying |FALSE|. The array therefore
% implies that at the beginning of the sequence, none of the points
% coincide with clasts. Later on in the |while| loop, we develop a new
% array called |ClastMapCheck|, which contains what ever was previously in
% |ClastMap|. We then apply the logical statement:
%
% |ClastMapCheck(r<rC)=ClastMapCheck(r<rC)+1;| 
%
% which adds one to all the element values in the |ClastMapCheck|, which are
% situated within the newly generated clast.
%
% We then perform a check, |max(ClastMapCheck(:))<=1|, which will be true
% if the new clast does not conicide with any existing clasts. If that is
% the case, the |ClastMap| array is updated with the values contained
% within |ClastMapCheck|, which means the new clast is incorporated into
% the collection. The clast counter, |ClastNo|, is then increased by one
% and so on.
%
% If the above check fails, this implies that the new clast coincides
% with existing clasts. The clast is therefore rejected and not
% incorporated into |ClastMap|.
%
% Now that we know which points, $(x_i,y_j)$, are located within the clasts
% we can add the following code to specify 2D arrays for the initial condition and the diffusion
% coefficient:
%
%   %Define initial condition
%   C0=zeros(size(x))+C0M;
%   C0(ClastMap==1)=C0C;
%   %Define diffusion coeffients
%   Da=zeros(size(x))+DaM;
%   Da(ClastMap==1)=DaC;
%
% Recall from the finite difference equations for the diffusive fluxes above, we
% require values at the block boundaries, |xB| and |yB|, notably,
% $D_{A,i+1/2,j}$ and $D_{A,i,j+1/2}$. Let us therefore introduce two new
% 2D arrays: |Dax| $=D_{A,i+1/2,j}$ and |Day| $=D_{A,i,j+1/2}$. These can
% be determined by interpolation using the code below:
%
%   %Interpolate the diffusive fluxes to the xB values
%   Dax=interp1(x(:,1),Da,xB(:,1));
%   %Interpolate the diffusive fluxes to the yB values
%   Day=interp1(y(1,:),Da',yB(1,:))';
% 
% Note the need to use the transpose operator, |'|, twice for the
% determination of |Day|. This is because |interp1| prefers to operate on columns. Another
% approach could be to use |interp2|. But using |interp1| like this could
% be more useful if we were to consider n-dimensional problems in the
% future.
%
% Next we need to add some code to specify the times that we want to solve for. We
% will ask for solutions at 0, 1, 10 and 100 years:
%
%   %Define times to solve for in years
%   tYR=[0 1 10 100];
%   %Convert times to seconds
%   t=tYR*365*24*60^2;
% 
% In the previous 1D PDE exercises, the resulting vector of ODEs (following
% the spatial discretisation) was solved using |ode45| and |ode15s|. However, 
% this 2D example leads to a 2D array of ODEs. MATLAB's ODE solvers 
% are only able to solve 1D vectors of ODEs. Therefore, the problem must be
% vectorised. This can be done using the |reshape| command by adding the
% following code:
%
%   %Vectorize the initial condition
%   C0=reshape(C0,Nx*Ny,1);
%
% This piece of code has transformed our 2D array, of $N_x$ rows and
% $N_y$ columns, into a 1D vector of $N_x\times N_y$ rows (and just one column).
% Read the help file for |reshape| to learn more.
%
% Now we will add some code to call |ode45| to solve the problem for us
%
%   %Ignore options for the ODE solver
%   options=[];
%
%   %Solve problem
%   h=waitbar(0,'Please wait...');
%   [t,C]=ode45(@MYodefun,t,C0,options,x,xB,y,yB,Nx,Ny,Dax,Day,h,t(end));
%   close(h)
%
% Because this solution will take a little longer to compile, we have also
% included an instruction to monitor the progress of the solver using a
% waitbar. The |h| term is a handle from the waitbar that is passed to the
% ODE function. It is also necessary to pass |t(end)| to remind the solver
% how long the simulation must run for. Read the help file for |waitbar| to
% learn more.
%
% Now we will write the ODE function. Copy the following code in a
% subfunction underneath your main function in the current script-file.
%
%   function dCdt=MYodefun(t,C,x,xB,y,yB,Nx,Ny,Dax,Day,h,tmax)
%
% The variable |C| contains a vector of all the $C_{i,j}$ values for the
% current time, |t|. Before any further work can be done with this, we need
% to de-vectorize this back to a 2D array containing $N_x$ rows and
% $N_y$ columns. This can be done using |reshape| by adding the
% following code.
%
%   %De-vectorize the C vector
%   C=reshape(C,Nx,Ny);
%
% Now we will calculate the diffusive fluxes for the internal block
% boundaries (i.e., not including the external boundaries).
%
%   %Calculate the diffusive fluxes
%   qx=-Dax(2:Nx,:).*diff(C,1,1)./diff(x,1,1);
%   qy=-Day(:,2:Ny).*diff(C,1,2)./diff(y,1,2);
%
% Two interesting things to note are: (1) We use a "2" in the |diff|
% input arguments for |qy| because we are now differencing the columns as
% opposed to the rows (study the help file for |diff| again if necessary);
% (2) we are only using |2:Nx| and |2:Ny| for |Dax| and |Day|,
% respectively, because we are not including the external boundaries at this stage.
%
% Now we will add some code to concatenate the zero-flux vectors to account for 
% the zero-flux boundary conditions.
%
%   %Apply zero flux boundary conditions
%   qx=[zeros(1,Ny); qx; zeros(1,Ny)];
%   qy=[zeros(Nx,1) qy zeros(Nx,1)];
%
% Next we can apply the mass conservation equation to calculate the
% $\partial C/\partial t|_{i,j}$ derivatives.
%
%   %Apply the mass conservation equation
%   dCdt=-diff(qx,1,1)./diff(xB,1,1)-diff(qy,1,2)./diff(yB,1,2);
%
% The above ODEs will be set in a 2D array, again of $N_x$ rows and
% $N_y$ columns. The following code is therefore needed to convert the 2D
% array back to a 1D vector before it can be sent back to the ODE solver.
%
%   %Vectorise the C vector
%   dCdt=reshape(dCdt,Nx*Ny,1);
%
% Finally, we will add the following code to monitor the progress of the
% solver using the waitbar.
%
%   %Monitor progress
%   waitbar(t/tmax,h)
%
%% Plotting the results as a set of 2D surfaces
%
% The following code can be added to plot the solution as a set of 2D
% surfaces.
%
%   %Plot the results
%   figure(1)
%   for n=1:4
%       subplot(2,2,n)
%       %Show results as a coloured surface
%       surf(x,y,reshape(C(n,:)',Nx,Ny))
%       %Set the x and y axis limits
%       axis([0 Lx 0 Ly])
%       %Set the color-scale limites
%       caxis([C0M C0C])
%       %Display the color bar
%       colorbar
%       %Use linear interpolation in the colour shading
%       shading interp
%       %Show a plan view of the surface
%       view(2)
%       %Label x and y axis
%       xlabel('Distance (cm)')
%       ylabel('Distance (cm)')
%       %Add a title to each subplot
%       title(['Iron oxide mass fraction after ' num2str(tYR(n)) ' years'])
%   end
%
% Note the need to reshape and transpose the solution vectors back to a 2D
% array in the |surf| call. Read the help file for |surf| if necessary.
%
%% Developing a more efficient solution using |ode15s|
%
% The problem should run reasonably fast. However, |ode45| is not the best
% solver for this problem. The different diffusion coefficients in
% the clasts and the matrix give rise to a stiff problem. Therefore, it
% would be better to use |ode15s|.
%
% But interestingly, you will find that simply
% changing the solver to |ode15s| makes the simulation take significantly
% longer. The reason for this is that |ode15s| is using finite differences
% to calculate the Jacobian matrix for an $N_x\times N_y$ system, which
% will be of size, $(N_xN_y)^2$. But the computation time can be significantly
% reduced by specifying the Jacobian pattern.
%
% Further consideration of the governing finite difference equations suggest 
%
% $\displaystyle \left.\frac{dc}{dt}\right|_{i}=a_{i,j}c_{i,j-1}+
% +b_{i,j}c_{i-1,j}+d_{i,j}c_{i,j}+e_{i,j}c_{i+1,j}+f_{i,j}c_{i,j+1}$
%
% from which we see that the Jacobian pattern has five non-zero diagonals.
%
% Modify the existing code (concerning |options=[]|) such that the Jacobian
% pattern is specified as follows:
%
%   %Define and set the Jacobian pattern
%   JPat=spdiags(ones(Nx*Ny,5),[-Nx -1 0 1 Nx],Nx*Ny,Nx*Ny);
%   options=odeset('JPattern',JPat);
%
% The |spdiags| command enables us to construct a five-diagonal sparse
% matrix. The |ones(Nx*Ny,5)| tells us there will be five diagonals containing
% ones. The |[-Nx -1 0 1 Nx]| prescribes the column number at which the diagonals
% will start in the top row of our sparse matrix. The |Nx*Ny,Nx*Ny| tells us that
% the resulting sparse matrix will be of size, $(N_xN_y)^2$.
%
% The code should run much faster now.