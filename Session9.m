%% NNESMO: Session 9 - Solving non-linear PDEs
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
% * Derive multiple finite difference solutions for the Burgers'
% equation.
% * Appreciate the difficulty of specifying stability criteria for non-linear problems.
% * Solve non-linear partial differential equations (PDE) using a semi-implicit scheme.
% * Use a Newton iteration method to develop a fully implicit method for non-linear PDEs.
% * Solve non-linear PDEs using MATLAB's ODE solver, ODE15s.
%
%% The Burgers' equation (BE)
%
% The Burgers' equation (BE) represents a  fluid momentum conservation statement
% in the presence of neglible pressure changes and body forces. 
% For a one-dimensional system, the BE can be written out as follows:
%
% $$ \displaystyle \frac{\partial u}{\partial t}+u\frac{\partial u}{\partial
% x}-\nu\frac{\partial^2 u}{\partial x^2}=0 $$ _____(1a)
%
% where $u$ [ $\mathrm{LT}^{-1}$ ] is velocity,  $x$ [ $\mathrm{L}$ ] is distance,
% $t$ [ $\mathrm{T}$ ] is time and $\nu$ [ $\mathrm{L}^2\mathrm{T}^{-1}$ ] is dynamic viscosity.
%
% The BE is interesting to study from a computational perspective because
% it represents a simple modification of the advection diffusion equation (ADE)
% such that a non-linear partial differential equation (PDE) is produced. 
% Comparing Eq. (1a) with Eq. (1) from the previous session it can be seen that  
% $a=u$ and $b=\nu$. Non-linear forms of the ADE come about in the Earth
% sciences for a variety of different applications. For example, fluid flow
% through unsaturated porous media (the Richards' equation), or heat
% conduction through rock formations where the material properties (e.g.,
% specific heat capacity and conductivity) have a strong dependence on
% temperature.
%
% As will be seen subsequently, adding just a small amount of non-linearity to
% a PDE leads to the need for quite different treatment when seeking to
% develop numerical solutions using approximate techniques such as finite
% differences.
%
% For the next few examples, the following initial and boundary conditions
% are considered:
%
% $$ \begin{array}{lll}
% u=u_0, & |x| \leq x_0, & t=0\\
% u=0, & x_0< |x| \leq L, & t=0\\
% u=0, & |x| = L, & t>0\\
% \end{array} $$ _____(1b)
%
% where $u_0$ [ $\mathrm{LT}^{-1}$ ] is the initial velocity of a region of
% $2x_0$ [ $\mathrm{L}$ ] length and $2L$ [ $\mathrm{L}$ ] represents the
% distance between the  two boundary conditions.
%
% For the special case where $L\rightarrow \infty$, the above system of
% equations can be solved exactly to obtain (McElwaine, 2016, pers. com.)
%
% $$ \displaystyle  \frac{u}{u_0}=\left[1+\frac{\displaystyle\exp\left(\frac{x+x_0-u_0t/2}{2\nu/u_0}\right)
%     \mathrm{erfc}\left(\frac{x_0+x}{2\sqrt{\nu t}}\right)
%     +\exp\left(\frac{x-x_0-u_0t/2}{2\nu/u_0}\right)
%     \mathrm{erfc}\left(\frac{x_0-x}{2\sqrt{\nu t}}\right)}
%     {\displaystyle\mathrm{erfc}\left(\frac{x-x_0-u_0t}{2\sqrt{\nu t}}\right)
%         -\mathrm{erfc}\left(\frac{x+x_0-u_0t}{2\sqrt{\nu
%         t}}\right)}\right]^{-1} $$  _____(1c)
%
% To evaluate the analytical solution, type in the following code in a single 
% script-file and save it in an appropriate place as "BurgersStudy.m".
%
%   function BurgersStudy
%   %A sequence of different subfunctions to solve the Burgers' equation.
%   
%   %Define model parameters
%   nu=0.1; %m^2/day
%   u0=1.2; %m/day
%   x0=0.7; %m
%
%   %Define location of points for plotting
%   x=linspace(-4,6,100);
%   %Define values of time (days) to be used for plotting
%   t=[0.5 1 2 4 8];
%   
%   %Evaluate the exact solution
%   uA=AnalSol(x,t,nu,u0,x0);
%   
%   %Plot results
%   figure(1)
%   hold off
%   plot(x,uA,'.')
%   xlabel('Distance (m)')
%   ylabel('Velocity (m/day)')
%   for n=1:numel(t)
%       MYlegend{n}=[num2str(t(n)) ' days'];
%   end
%   legend(MYlegend,'location','northwest')
%   
%   %**********************************************************************
%   
%   function u=AnalSol(x,t,nu,u0,x0)
%   %Subfunction containing the analytical solution
%   [x,t]=ndgrid(x,t);
%   x1=(x+x0-u0*t)./sqrt(4*nu*t);
%   x2=(x-x0-u0*t)./sqrt(4*nu*t);
%   x3=(x0+x)./sqrt(4*nu*t);
%   x4=(x0-x)./sqrt(4*nu*t);
%   c=u0/2/nu;
%   T3=exp(c*(x+x0-u0*t/2));
%   T4=exp(c*(x-x0-u0*t/2));
%   u=u0./(1+(T3.*erfc(x3)+T4.*erfc(x4))./(erfc(x2)-erfc(x1)));
%
% The proceeding outline of this session is as follows. First 
% a simple explicit time-stepping scheme for the Burgers' equation is
% devloped and  its stability requirement is assessed. Following on from this, a more
% stable implicit scheme is developed. However, because of the
% non-linearity in the Burgers' equation, the implicit scheme cannot be
% implemented directly. Therefore, the Newton iteration method is
% introduced as a methodology for iteratively solving the implicit scheme.
% Subsequently it is shown that the first iteration of the fully implicit
% scheme gives rise to an alternative unconditionally stable scheme,
% referred to as a semi-implicit scheme. Finally, there is a classroom
% assignment, which involves modifying code from "ADEStudy.m" (developed
% during the previous session) to evaluate the various possible finite
% difference schemes that can be used to solve the Burgers' equation.
%
%% ETCS - Explicit time-stepping with central differences in space
%
% First we will develop an Euler explicit time-stepping scheme for the BE 
% problem described above and use a
% central difference approximation for the $\partial u/\partial x$ term.
% Ignoring all the truncation errors associated with the Taylor expansions
% (discussed in the previous session),
% the BE can then be written as
%
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t} +u_{i,n}\left(\frac{u_{i+1,n}-u_{i-1,n}}
% {2 \Delta x}\right)
% -\nu\left(\frac{u_{i+1,n}-2u_{i,n}+u_{i-1,n}}{\Delta x^2}\right)=0 $$ 
%
% which can be rearranged to get
%
% $$ u_{i,n+1}=A_{i,n} u_{i-1,n}+B_{i,n} u_{i,n} + C_{i,n} u_{i+1,n} $$
% 
% where
% 
% $$ \displaystyle A_{i,n}=\frac{u_{i,n}\Delta t}{2 \Delta x}+\frac{\nu\Delta t}{\Delta x^2},
% \quad
% B_{i,n}=1-\frac{2\nu\Delta t}{\Delta x^2},
% \quad
% C_{i,n}=-\frac{u_{i,n}\Delta t}{2 \Delta x}+\frac{\nu\Delta t}{\Delta x^2} $$
%
% or alternatively it can be said that
%
% $$ u_{i,n+1}=P_{i,n} u_{i-1,n}+Q_{i,n} u_{i,n} + R_{i,n} u_{i+1,n} $$
% 
% where
%
% $$ \displaystyle P_{i,n}=\frac{\nu\Delta t}{\Delta x^2},
% \quad
% Q_{i,n}=1-\frac{2\nu\Delta t}{\Delta x^2}-\Delta t \left(\frac{u_{i+1,n}-u_{i-1,n}}{2 \Delta x}\right),
% \quad
% R_{i,n}=\frac{\nu\Delta t}{\Delta x^2} $$
% 
% Hereafter, both of the above schemes will be referred to as the ETCS scheme.
%
% Recall that for numerical stability of explicit time-stepping schemes to be ensured, it is
% generally found that the coefficients, $A_{i,n}$, $B_{i,n}$ and $C_{i,n}$ (or $P_{i,n}$, $Q_{i,n}$ and $R_{i,n}$), need
% always to be positive. For the case when $u$, $\partial u/\partial x$, $\nu$, $\Delta x$ and $\Delta
% t$ are positive, it follows that there are two
% stability criteria for these schemes:
% 
% $$ \displaystyle \mathrm{Cr_{eff}}\equiv \frac{\nu\Delta t}{\Delta x^2} +\frac{\Delta t}{2} \left(\frac{u_{i+1,n}-u_{i-1,n}}{2 \Delta x}\right)<0.5 
% \quad
% \mathrm{and}
% \quad
% \mathrm{Pe_{eff}}\equiv  \frac{u_{i,n}\Delta x}{\nu }<2 $$
% 
% where $\mathrm{Cr_{eff}}$ and $\mathrm{Pe_{eff}}$ can be thought of as being
% an effective Courant number and effective Peclet number, respectively.
%
% The issue here, which is common to most non-linear PDEs, is
% that it is not straightforward to specify a priori a $\Delta x$ and $\Delta t$
% that will ensure stability of the system for the desired simulation
% duration.
%
% Considering our experience from the previous session, it is postulated
% that an impliciting time-stepping scheme with a backward difference for
% the $\partial u/\partial x$ is likely to lead to an unconditionally
% stable scheme.
%
%% ITBS - Implicit time-stepping with backward differences
%
% If we consider an Euler implicit time-stepping scheme and use a
% backward difference approximation for the $\partial u/\partial x$ term,
% ignoring all the truncation errors associated with the Taylor expansions,
% the BE can be written as
%
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t}
% +u_{i,n+1}\left(\frac{u_{i,n+1}-u_{i-1,n+1}}{ \Delta
% x}\right)-\nu\left(\frac{u_{i+1,n+1}-2u_{i,n+1}+u_{i-1,n+1}}{\Delta
% x^2}\right)=0 $$
%
% which can be rearranged to get
%
% $$ u_{i,n}=A_{i,n+1} u_{i-1,n+1}+B_{i,n+1} u_{i,n+1} + C_{i,n+1} u_{i+1,n+1} $$ _____(3a)
% 
% where
% 
% $$ \displaystyle A_{i,n+1}=-\frac{u_{i,n+1}\Delta t}{ \Delta x}-\frac{\nu\Delta t}{\Delta x^2},
% \quad
% B_{i,n+1}=1+\frac{u_{i,n+1}\Delta t}{ \Delta x}+\frac{2\nu\Delta t}{\Delta x^2},
% \quad
% C_{i,n+1}=-\frac{\nu\Delta t}{\Delta x^2} $$
%
% or alternatively it can be said that
%
% $$ u_{i,n}=P_{i,n+1} u_{i-1,n+1}+Q_{i,n+1} u_{i,n+1} + R_{i,n+1} u_{i+1,n+1} $$
%
% where
% 
% $$ \displaystyle P_{i,n+1}=-\frac{\nu\Delta t}{\Delta x^2},
% \quad
% Q_{i,n+1}=1+\Delta t\left(\frac{u_{i,n+1}-u_{i-1,n+1}}{ \Delta x}\right)+\frac{2\nu\Delta t}{\Delta x^2},
% \quad
% R_{i,n+1}=-\frac{\nu\Delta t}{\Delta x^2} $$
%
% Hereafter this scheme will be referred to as the ITBS scheme.
%
% Recalling from the previous session, for implicit time-stepping
% schemes, the coefficients $A_{i,n+1}$ and $C_{i,n+1}$ (or $P_{i,n+1}$ and $R_{i,n+1}$) 
% need to be negative and the coefficents $B_{i,n+1}$ (or $Q_{i,n+1}$)
% need to be positive to ensure stability. Consequently, for the case when $u$, 
% $\partial u/\partial x$, $\nu$, $\Delta x$ and $\Delta t$ are positive,
% the ITBS scheme should be unconditionally stable. However, for the
% situation dicated by the boundary conditions given in Eq. (1b), even for
% positive $u_0$, there will be a region where $\partial u/\partial x<0$.
%
% Despite the fact that ITBS does not appear to be unconditionally stable,
% it can be appreciated that this scheme is likely to be much more stable
% than the ETCS scheme.
%
% However, recall that the $\mathbf{A}$, $\mathbf{B}$ and $\mathbf{C}$ vectors of coefficients
% represent the first, second and third diagonals of a tridiagonal Jacobian
% matrix, $\mathbf{J}_{n+1}$, defined by the formula
%
% $$ \mathbf{u}_n=\mathbf{J}_{n+1} \mathbf{u}_{n+1} $$ _____(3b)
%
% Another challenge here is that evaluation of $\mathbf{J}_{n+1}$
% also requires values from the $\mathbf{u}_{n+1}$ vector. Consequently
% Eq. (3b) must be solved iteratively.
%
%% The Newton iteration method
%
% The Newton iteration method (also referred to as the Newton-Raphson
% method), is a simple root finding algorithm often used to solve
% non-linear PDEs using implicit time-stepping schemes such as Eq. (3b) above.
%
% The Newton iteration method for root finding can be simply derived as
% follows:
%
% Consider a function, $f(x)$. It is desired to find a value, $x_{\mathrm{root}}$, such that
% $f(x_{\mathrm{root}})=0$ (i.e., a root of $f(x)$) that is situated close to a value $x_0$.
%
% The equation for a line, $y$, that lies tangent to $f(x)$ at $x=x_0$ can
% be written as
%
% $$ y = f'(x_0) (x-x_0) +f(x_0) $$ _____(4a)
%
% A good estimate for a value of $x_{root}$, which we will call $x_1$, can be found by extrapolating
% the tangent equation, Eq. (4a), to $y=0$, i.e.:
%
% $$ 0 = f'(x_0) (x_1-x_0) +f(x_0) $$
%
% which can be rearranged to get
%
% $$ x_1 = x_0 - f(x_0) / f'(x_0) $$
%
% from which it can be envisaged that a progressively closer estimate of $x_{root}$
% can be obtained from the recursive relationship
%
% $$ x_{m+1} = x_m - f(x_m) / f'(x_m) $$ _____(4b)
%
% where $m$ denotes the number of iterations that have been undertaken.
%
% The above system can also be used to solve systems of non-linear equations
% such as Eq. (3b). In this case it is desired to find the roots of
% the vector, $\mathbf{F}_{n+1,m}$, where $m$ denotes results from the m-th 
% iteration, which satisfies
%
% $$  \mathbf{F}_{n+1,m} = \mathbf{u}_n-\mathbf{J}_{n+1,m} \mathbf{u}_{n+1,m} $$
%
% from which it can be seen that 
%
% $$\displaystyle \frac{\partial \mathbf{F}_{n+1,m}}{\partial \mathbf{u}_{n+1,m}} =-\mathbf{J}_{n+1,m}$$
%
% Further consideration of Eq. (4b) then suggests that a progressively closer estimate of $\mathbf{u}_{n+1}$
% can be obtained from 
%
% $$ \mathbf{u}_{n+1,m+1} = \mathbf{u}_{n+1,m} +\mathbf{J}_{n+1,m}^{-1} [\mathbf{u}_n-\mathbf{J}_{n+1,m} \mathbf{u}_{n+1,m}] $$
%
% which reduces to
%
% $$ \mathbf{u}_{n+1,m+1} = \mathbf{J}_{n+1,m}^{-1} \mathbf{u}_n $$ _____(4c)
%
% where $\mathbf{J}_{n+1,m}$ is the Jacobian matrix (as defined in Eq.
% (3b)) as determined assuming $\mathbf{u}_{n+1}=\mathbf{u}_{n+1,m}$.
%
% A good starting point is to try $\mathbf{u}_{n+1,0}=\mathbf{u}_{n}$.
%
% The iterative process should proceed until the scheme has achieved an
% acceptable level of local convergence. A suitable criterion for the
% accepting of $\mathbf{u}_{n+1,m}$ as a good enough approximation of $\mathbf{u}_{n+1}$
% can be written as follows:
%
% $$\displaystyle \left|1 - \frac{u_{i,n+1,m}}{u_{i,n+1,m+1}}\right| < 10^{-6}, \quad i=1\ldots N_x$$
%
% where $N_x$ denotes the number of points that $x$ has been discretised
% into.
%
%% SITBS - Semi-implicit time-stepping with backward differences
%
% Note that if only one iteration of the above scheme is considered we have
%
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t}
% +u_{i,n}\left(\frac{u_{i,n+1}-u_{i-1,n+1}}{ \Delta
% x}\right)-\nu\left(\frac{u_{i+1,n+1}-2u_{i,n+1}+u_{i-1,n+1}}{\Delta
% x^2}\right)=0 $$
%
% which can be rearranged to get
%
% $$ u_{i,n}=A_{i,n+1} u_{i-1,n+1}+B_{i,n+1} u_{i,n+1} + C_{i,n+1} u_{i+1,n+1} $$ _____(5a)
% 
% where
% 
% $$ \displaystyle A_{i,n+1}=-\frac{u_{i,n}\Delta t}{ \Delta x}-\frac{\nu\Delta t}{\Delta x^2},
% \quad
% B_{i,n+1}=1+\frac{u_{i,n}\Delta t}{ \Delta x}+\frac{2\nu\Delta t}{\Delta x^2},
% \quad
% C_{i,n+1}=-\frac{\nu\Delta t}{\Delta x^2} $$
%
% which, providing $u_{i,n}>0$, can be seen to be unconditionally stable.
%
% This scheme should be hereafter referred to as SITBS.
%
% Such a scheme is often referred to as a semi-implicit scheme or a mixed
% implicit and explicit scheme because it is implicit in the sense that
% Eq. (5a) cannot be solved explicitly for $u_{i,n+1}$, but explicit in
% terms of the fact that the $u_i$ terms in the $A_{i,n}$ and
% $B_{i,n}$ coefficients are taken from the previous time-step. Such
% schemes are often implemented as a simple alternative to  fully implicit
% schemes in practice.
%
%% Classroom Assignment
%
% 1) Modify the code developed in ADEStudy.m (developed during the previous
% session) to add an appropriate subfunction to BurgersStudy.m to 
% solve the problem described by Eqs. (1a) and (1b) using
% the SITBS scheme and compare your results to the analytical solution
% given in Eq. (1c). Try using $\Delta x=0.02$ m and $\Delta t=0.25$ days.
% Take special care when incorporating the $\mathbf{A}$, $\mathbf{B}$ and
% $\mathbf{C}$ vectors into your Jacobian matrix using |spdiags|. Check the 
% help documentation for |spdiags| to be sure.
%
% 2) Modify your SITBS scheme further to obtain the ITBS scheme using the
% Newton iteration method described above. Specify a maximum number of
% iterations of 100. How many iterations are required for convergence? 
% How much more accurate is ITBS as compared SITBS?
%
% 3) Solve the same problem using |ODE15s| with backward differencing. Note that specifying the
% Jacobian pattern (as we did in previous sessions) should lead to
% a considerable CPU time saving. How does this solution compare to ITBS?
% What happens if you set $\nu=0.01$ $\mathrm{m}^2\mathrm{day}^{-1}$ ? Try
% using |tic| and |toc| (see help documentation) to determine which of your 
% approximate solutions is most
% efficient? What happens if you use a central differencing scheme instead?

    

