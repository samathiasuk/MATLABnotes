%% MATLAB Notes: Session 6 - Solving PDEs using ODE solvers
%
% Simon Mathias
%
% Department of Engineering
%
% Durham University
%
% <https://samathiasuk.github.io/MATLABnotes/html/MATLABcontents.html Return to main contents page>
%
%% Learning outcomes
%
% At the end of the session you should be able to:
%
% * Use a finite difference spatial discretisation to transform a partial
% differential equation (PDE) into a set of coupled ordinary differential equations (ODE).
% * Solve the one-dimensional advection dispersion equation using |ode45| or |ode15s|.
% * Apply a non-uniform grid spacing.
% * Understand the meaning of stiff and non-stiff problems.
% * Determine the Jacobian pattern associated with a PDE problem.
% * Solve PDEs using |ode15s|.
%
%% Groundwater flow in a confined aquifer
%
% A confined aquifer is bounded to the West by a lake and to the East by an 
% impermeable fault zone. The aquifer is homogenous and isotropic and 
% characterised by a transmissivity of 800 $\mathrm{m}^{2}\mathrm{day}^{-1}$
% and a storativity of 0.01. The edge of the lake lies parallel to the 
% fault zone and is separated by a distance of 1600 m. Initially the lake 
% level is 50  mAOD (metres above ordinance datum). After a significant episode of snow melt in the 
% mountains above, the water level in the lake is suddenly raised to 53
% mAOD.
% 
% The governing equation for one dimensional groundwater flow in a confined
% aquifer takes the form
%
% $$ \displaystyle S\frac{\partial h}{\partial t}=-\frac{\partial Q}{\partial x} $$
% 
% where $S$ [-] is the aquifer storativity (which relates to the bulk
% compressibility of the aquifer), $h$ [L] is the hydraulic head, $t$ [T]
% is time, $x$ [L] is distance and $Q$ $[\mathrm{L}^2\mathrm{T}^{-1}]$
% is the volumetric flow rate per unit breadth of confined aquifer, found
% from
%
% $$ \displaystyle Q=-T\frac{\partial h}{\partial x} $$
%
% where $T$ $[\mathrm{L}^2\mathrm{T}^{-1}]$ is the transimssivity of
% the aquifer (which relates to the bulk permeability of the aquifer).
%
% The relevant initial and boundary conditions can be written as follows:
%
% $$ \begin{array}{lll}
% h=h_I, & 0\leq x \leq L, & t=0\\
% h=h_0, & x=0, & t>0\\
% Q=0, & x=L, & t>0\\
% \end{array} $$
%
% where in this case, $h_I=50$ m, $h_0=53$ m and $L=1600$ m.
%
% Use |ode45| to develop a one-dimensional finite difference model to estimate 
% the hydraulic head within the aquifer.
% 
% Show your results as a plot of hydraulic head against distance, $x$, for 
% the following times: 0.1, 0.3, 0.5, 1, 3, 5 and 7 days.
%
%% Spatial discretisation using finite differences
%
% The above problem is an example of a partial differential equation (PDE).
% However, if we discretise the problem in space, the problem becomes a
% coupled set of ordinary differential equations (ODE) with respect to
% time. In the previous exercise we used |ode45| to solve a single ODE.
% Here we will use |ode45| to solve the large set of coupled ODEs that
% results from the spatial discretisation. We will use finite difference 
% for spatial discretisation. However, other methods such as finite
% elements and pseudospectral methods can also be used in a similar way.
%
% Let us consider a set of $N$ discrete points on the $x$-axis:
% $x_1,x_2,x_3,\ldots,x_{N-1}, x_N$.
%
% The corresponding set of hydraulic heads can be written as: $h_1,h_2,h_3,\ldots,h_{N-1}, h_N$.
%
% An approximation of flow per unit breadth, $Q$, can be found from
%
% $$\displaystyle Q_{i+1/2}=-T\left(\frac{h_{i+1}-h_i}{x_{i+1}-x_i}\right) $$
%
% but note that these flow rates are defined at an alternative set of
% points:  $x_{1/2},x_{1+1/2},x_{2+1/2},\ldots,x_{N-1/2}, x_{N+1/2}$. This
% was discussed to some extent in "NNESMO: Session 4 -
% Approximate methods for differentiation and integration".
%
% The resulting set of ODEs to be solved takes the form
%
% $$ \displaystyle \left.\frac{dh}{dt}\right|_{i}=-\frac{1}{S}\left(
% \frac{Q_{i+1/2}-Q_{i-1/2}}{x_{i+1/2}-x_{i-1/2}}\right) $$
%
% Recall that the relationship between $x_{i+1/2}$ and $x_i$ is
%
% $$ x_{i+1/2}=(x_{i+1}+x_i)/2 $$
%
% Furthermore it can be shown that
%
% $$ x_{i}=(x_{i+1/2}+x_{i-1/2})/2 $$
%
%% Development of the solution code
%
% Create a new script-file and type the following:
%
%   function MATLABsession6_Assignment
%   %This file contains instructions requested by MATLAB Notes - Session 6
%
%   %Compute instructions associated with the groundwater flow example
%   GroundwaterExample
%
%   %**************************************************************************
%
%   function GroundwaterExample
%   %Solves the one-dimensional confined aquifer flow equation using a MATLAB
%   %ODE solver
%
%   %T (m^2/day) - Transmissivity
%   %S (-) - Storativity
%   %L (m) - Length of the aquifer
%   %hI (m) - Initial hydraulic head
%   %h0 (m) - Boundary head
%   %x (m) - location of points being solved for
%   %xB (m) - location of finite difference block boundaries
%   %t (days) - time
%   
%   %Define model parameters
%   T=800;
%   S=0.01;
%   L=1600;
%   hI=50;
%   h0=53;
%
% So far all we have done is written some comments explaining what some of
% the subsequent variables are and defined the model parameter values.
%
% The next step is to determine the locations on the x-axis where we are
% going to solve for.
%
% In the above variable list we have |x| and |xB|. It is planned that |x| 
% and |xB| will be  a vectors containing 
% $x_1,x_2,x_3,\ldots,x_{N-1}, x_N$ and 
% $x_{1/2},x_{1+1/2},x_{2+1/2},\ldots,x_{N-1/2}, x_{N+1/2}$, respectively.
%
% Another point to note is that |xB(1)=0| and |xB(N+1)=L|.
%
% Add the following code to provide a discretisation for 20 equally spaced
% solution points  
%
%   %Define spatial grid
%   N=20;
%   xB=linspace(0,L,N+1)';
%   x=(xB(1:N,1)+xB(2:N+1,1))/2;
%
% Next we will add some code to define the times of interest and a vector
% containing the initial values for $h_i$.
%
%   %Define the times of interest
%   t=[0 0.1 1 10 100];
%
%   %Define intial condition vector for the ode solver
%   hI_vec=zeros(size(x))+hI;
%
%
% Now we will call |ode45| to provide our solution. Note that we also need
% to pass all the relevant input parameters to the ODE function, |MYodefun|.
%
%   %Apply ode45 to obtain the finite difference solution
%   options=[];
%   [t,hFD]=ode45(@GroundwaterODEfun,t,hI_vec,options,x,xB,S,T,h0);
%
%
% Of course we need to write the ODE function as well. This can be added as
% a subfunction underneath the |GroundwaterExample| function.
%
%   function dhdt=GroundwaterODEfun(t,h,x,xB,S,T,h0)
%   %Calculate derivatives and include h=h0 at x=0 boundary condition
%   dhdx=diff([h0;h],1,1)./diff([xB(1);x],1,1);
%   %Calculate Darcy fluxes and include Q=0 at x=L boundary condition
%   Q=[-T*dhdx;0];
%   %Calculate the flux divergence
%   dQdx=diff(Q,1,1)./diff(xB,1,1);
%   %Solve the mass conservation statement
%   dhdt=-dQdx/S;
%
%
% Most of the above code in the ODE function is self-explanatory. However, 
% pay close attention
% to how the boundary conditions are applied.
%
% The fixed head boundary, associated with the lake,
% is applied by concatenating the boundary head, |h0|, and its associated
% location, |xB(1)|, to the |h| and |x| vectors, respectively, prior to
% finding the derivatives, $dh/dx$.
%
% The zero flux boundary, associated with the fault zone, is applied by
% simply concatenating a zero to the |Q| vector.
%
%% Plotting the results and comparing to an analytical solution
%
% When ever developing numerical models it is always important to study your
% results and try to compare these to analytical solutions where possible.
%
% The problem being solved here has an analytical solution for the special
% case where $L\rightarrow\infty$. Our finite difference model should
% produce very similar results to the analytical solution until the
% pressure wave hits the bounday at $x=L$.
%
% The analytical solution takes the form
%
% $$ \displaystyle h=(h_0-h_I)\mathrm{erfc}\left(\sqrt{\frac{Sx^2}{4Tt}}\right)+h_I $$
%
% Add the following subfunction to your script-file, which contains an
% implementation of the analytical solution.
%
%   function h=GroundwaterAnalyticalSol(t,x,S,T,hI,h0)
%   [t,x]=ndgrid(t,x);
%   h=(h0-hI)*erfc(sqrt(S*x.^2/4/T./t))+hI;
%
% Add the following code to the |GroundwaterExample| function in your script-file to
% evaluate the analytical solution and generate a plot comparing the results 
% from both the analytical and finite difference solutions.
%
%   %Evaluate analytical solution
%   h=GroundwaterAnalyticalSol(t,x,S,T,hI,h0);
%   
%   %Plot results and compare with the analytical solution
%   figure(1)
%   clf
%   hold on
%   plot(x,h)
%   %Reset colour order for plotting
%   set(gca,'ColorOrderIndex',1)
%   plot(x,hFD,'o--')
%   xlabel('Distance (m)')
%   ylabel('Hydraulic head (m)')
%   legend('0 days (analytical)','0.1 days (analytical)',...
%       '1 days (analytical)','10 days (analytical)','100 days (analytical)'...
%       ,'0 days (finite difference)','0.1 days (finite difference)'...
%       ,'1 days (finite difference)','10 days (finite difference)'...
%       ,'100 days (finite difference)')
%   %Make a nice box around the graphs
%   box on
%
%% Stiff and non-stiff problems
%
% At the moment, we are solving the diffusion problem using a uniform
% space-step. Just as we can benefit from having time-step refinement
% during periods of high activity, we can also benefit from spatial grid
% refinement in areas of high activity.
%
% In this particular case, most of the activity is occuring around $x=0$ where
% the fixed head boundary exists. Therefore a more appropriate spatial
% discretisation is arguably achieved using the code:
%
%   xB=[0;logspace(-3,0,N)'*L];
%
% which logarithmically spaces the finite different points over three
% orders of magnitude, with the smallest grid spacing around $x=0$. Look up
% |logspace| in the help file to find out more.
%
% Add the above code for |xB| to the existing code and run.
%
% You will find the solution takes forever to complete. The reason is that
% the simulation has become very stiff. Type "ctrl c" in the command
% window to terminate the simulation.
%
% Sets of coupled ODEs are said to be 
% stiff when the ODEs move at very different rates. By refining the spatial
% grid, we have created a system whereby the ODEs near $x=0$ are very fast
% whilst the ODEs near $x=L$ are very slow.
%
% The |ode45| solver is not good for stiff problems. Instead, try and use
% the solver |ode15s| by changing the code to say
%
%   [t,hFD]=ode15s(@GroundwaterODEfun,t,hI_vec,options,x,xB,S,T,h0);
%
% You can read more about how |ode15s| works in the help files.
%
%% Specifying the Jacobian pattern
%
% For large sets of ODEs, the stiff solvers are much more efficient if you
% specify the Jacobian pattern a priori. So what is the Jacobian pattern?
%
% For the example under consideration, let us define
%
% $$ \displaystyle f_i=\left. \frac{\partial h}{\partial t}\right|_{i} $$
%
% Now consider the vectors $\mathbf{f}=[f_1,f_2,f_3,\ldots,f_{N-1},f_{N}]^T$ and
% $\mathbf{h}=[h_1,h_2,h_3,\ldots,h_{N-1},h_{N}]^T$.
%
% There is a matrix, $\mathbf J$, that exists such that
%
% $$ \mathbf{f}=\mathbf{J} \mathbf{h} $$
%
% After some further consideration, it can be understood that this matrix
% is defined as follows:
%
% $$ \displaystyle \mathbf{J}=\frac{d\mathbf{f}}{d\mathbf{h}}=\left[
% \begin{array}{ccc}\displaystyle
% \frac{\partial f_1}{\partial h_1} & \ldots & \displaystyle\frac{\partial f_1}{\partial h_N}\\
% \vdots & \ddots & \vdots\\
% \displaystyle\frac{\partial f_N}{\partial h_1} & \ldots & \displaystyle\frac{\partial f_N}{\partial h_N}\\
% \end{array}
% \right] $$
%
% The $\mathbf J$ matrix quantatively described how each of the $f_i$
% values depend on each of the $h_i$ values. Such a type of matrix is often
% called a Jacobian matrix.
%
% By default, the stiff solvers in MATLAB calculate the Jacobian matrix 
% using a set of finite difference calculations. For a set of $N$ ODEs, the
% Jacobian matrix will have $N^2$ elements. Therefore, it can be understood
% that such an approach can become computationally very expensive. However,
% in practice, the Jacobian matrix is very sparse. Consequently far fewer
% calculations are necessary. Therefore, there are great benefits to be had
% from informing MATLAB where the location of non-zero values are. This
% can be done by specifying the so-called Jacobian pattern.
%
% Consider again the finite difference approximations for $Q$ and $\partial
% h/\partial t$ above. Substituting the equation for $Q_{i+1/2}$ into the
% equation for $dh/dt|_{i}$ leads to
%
% $$\displaystyle \left.\frac{dh}{dt}\right|_{i}=a_i h_{i-1}+b_ih_i+c_ih_{i+1} $$
%
% where
%
% $$ \displaystyle a_i = \frac{T}{S(x_{i}-x_{i-1})(x_{i+1/2}-x_{i-1/2})} $$
%
% $$\displaystyle c_i = \frac{T}{S(x_{i+1}-x_i)(x_{i+1/2}-x_{i-1/2})} $$
%
% $$\displaystyle b_i = -(a_i + c_i) $$
%
% Therefore, the Jacobian pattern for this problem can be seen to be a
% tri-diagonal sparse matrix of the form:
%
% $$ \left[
% \begin{array}{ccccccccc}
% 1 & 1 & 0 & 0 & \ldots & 0 & 0 & 0 & 0\\ 
% 1 & 1 & 1 & 0 & \ldots & 0 & 0 & 0 & 0\\ 
% 0 & 1 & 1 & 1 & \ldots & 0 & 0 & 0 & 0\\ 
% 0 & 0 & 1 & 1 & \ldots & 0 & 0 & 0 & 0\\ 
% \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots & \vdots\\ 
% 0 & 0 & 0 & 0 & \ldots & 1 & 1 & 0 & 0\\ 
% 0 & 0 & 0 & 0 & \ldots & 1 & 1 & 1 & 0\\ 
% 0 & 0 & 0 & 0 & \ldots & 0 & 1 & 1 & 1\\ 
% 0 & 0 & 0 & 0 & \ldots & 0 & 0 & 1 & 1\\ 
% \end{array}
% \right] $$
%
% To implement this Jacobian pattern within your code, replace the line
%
%   options=[];
%
% with
%
%   %Define and set the Jacobian pattern
%   JPat=spdiags(ones(N,3),[-1 0 1],N,N);
%   options=odeset('JPattern',JPat);
%
% Run your programme again. It should be a little faster. Specifying the
% Jacobian pattern will much more important when we look at 
% a two-dimensional problem in a later session.
%
% To understand further how, the Jacobian pattern has been specified, read
% the help files for |odeset| and then for |spdiags|.
%
% The |spdiags| command enables us to construct a tri-diagonal sparse
% matrix. The |ones(N,3)| tells us there will be three diagonals containing
% ones. The |[-1 0 1]| presribes the column number at which the diagonals
% will start in the top row of our sparse matrix. The |N,N| tells us that
% the resulting sparse matrix will be of size, $N\times N$.
%
%% Classroom assignment
%
% An incident at a paper factory has led to the dumping of product in a 
% nearby river. Creat a new subfunction called |ChemicalTransportExample| and use the advection dispersion equation to simulate the 
% subsequent migration of contamination.
%
% The relevant governing equation is:
%
% $$ \displaystyle \frac{\partial C}{\partial t}=-\frac{\partial q}{\partial x}$$
%
% where
%
% $$ \displaystyle q=vC-D_L\frac{\partial C}{\partial x} $$
%
% and $C$ $[\mathrm{ML}^{-3}]$ is the solute concentration of the
% product, $t$ $[\mathrm{T}]$ is time, $x$ [L] is space and $D_L$
% $[\mathrm{L}^2\mathrm{T}^{-1}]$ is the longitudinal dispersion coefficent.
%
% The relevant initial and boundary conditions can be written as follows:
%
% $$ \begin{array}{lll}
% C=C_I, & 0\leq x \leq L, & t=0\\
% C=C_0(t), & x=0, & t>0\\
% \displaystyle \partial C/\partial x=0, & x=L, & t>0\\
% \end{array} $$
%
% where in this case, $L=10$ km, $C_I=0.1$ mg/l and
%
% $$\displaystyle C_0(t)=\left\{\begin{array}{ll}
% C_{00}, & 0\leq t\leq t_{0}\\
% \\
% C_I, & t>t_{0}
% \end{array}\right.$$
%
% where $C_{00}=120$ mg/l and $t_{0}=1$ day.
%
% The velocity in the river is 1.2 km/day and the associated longitudinal
% dispersion coefficient is 0.23 $\mathrm{km}^2/\mathrm{day}$.
%
% Discretise your domain such that you have 50 equally-spaced solution
% points. Then write some code to solve the above problem using |ode15s|.
%
% Note that something interesting has to be done to obtain an appropiate
% vector of concentration values, which can be multplied by $v$ in the
% expression for $q$ above. Let us call this vector of concentration values
% |CB| (i.e., concentrations at the block boundaries). I suggest you adopt
% the following approach:
%
%   %Interpolate C to the xB points
%   CB=interp1(x,C,xB);
%   %Impose fixed concentration boundary
%   CB(1)=C0;
%   %Impose zero gradient boundary
%   CB(end)=C(end);
%
% Such an approach essentially represents a central differencing scheme.
% Read the help files about |interp1| if you have not seen this before.
% This is a very useful function.
%
% Paste the following subfunctions at the end of your code. These contain
% an analytical solution to the above problem:
%
%   function C=ChemicalTransportAnalyticalSol(t,x,v,DL,CI,C00,t0)
%   %Analytical solution for the advection dispersion problem
%   %t (days) - time after start of spill
%   %x (km) - distance from spill
%   %v (km/day) - velocity of river flow
%   %DL (km^2/day) - longitudinal dispersion coefficient
%   %CI (mg/l) - initial concenctration
%   %C00 (mg/l) - concentration at spill site during spill incident
%   %t0 (days) - duration of spill incident
%   [t,x]=ndgrid(t,x);
%   F1=FFun(t*v./x,v*x/DL);
%   F2=FFun((t-t0)*v./x,v*x/DL);
%   F2(t<t0)=0;
%   F=F1-F2;
%   C=(C00-CI)*F+CI;
%
%   %**************************************************************************
%   function F=FFun(z,P)
%   z(z<0)=0;
%   F=(erfc(sqrt(P/4./z).*(1-z))+exp(P).*erfc(sqrt(P/4./z).*(1+z)))/2;
%   ind=P>700;
%   F(ind)=erfc(sqrt(P(ind)/4./z(ind)).*(1-z(ind)))/2;
%
% Generate a plot of solute concentration against distance for 0, 2, 4 and
% 6 days. Compare the results from your finite difference code with those
% from the analytical solution.
%
% You should find there is almost a perfect correspondence between the
% finite difference code and the analytical solution. Now reduce the
% dispersion coefficient by a factor of 20. What happens and why?
%
% An example MATLAB code, containing all of the instructions requested
% above, is given in <https://samathiasuk.github.io/MATLABnotes/MATLABsession6_Assignment.m MATLABsession6_Assignment.m>.
 
