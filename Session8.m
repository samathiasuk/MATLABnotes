%% NNESMO: Session 8 - Stability and numerical diffusion
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
% * Derive multiple finite difference solutions for the advection diffusion
% equation.
% * Speculate about stability criteria by assessing the sign of finite difference 
% coefficients.
% * Calculate numerical diffusion due to different finite difference schemes
% using Taylor expansions.
%
%% The advection diffusion equation (ADE)
%
% Most of the problems considered up till now can be thought of as
% representing special forms of the so-called advection diffusion equation
% (ADE)
%
% $$ \displaystyle \frac{\partial u}{\partial t}+a\frac{\partial u}{\partial
% x}-b\frac{\partial^2 u}{\partial x^2}=0 $$ _____(1)
%
% where $u$ is the dependent variable, $x$ is distance,
% $t$ is time and $a$ and $b$ are constant parameters.
%
% We specifcally solved the ADE when
% considering the previous pollution transport scenario. Recall it was
% found that the ADE exhibited some numerical instabilities when high
% values of $a$ were applied. The reason that the ADE is challenging to
% solve in this context  
% is due to its composite hyperbolic and parabolic nature.
%
% When $a=0$, the ADE reduces to a parabolic partial differential equation
% (PDE) identical to the groundwater flow equation, chemical diffusion 
% equation and heat conduction equation considered previously. A parabolic
% equation is characterized by having all time derivatives of first-order
% and all spatial derivatives of second-order.
% Such equations are relatively
% straightforward to solve using finite difference methods.
%
% When $b=0$, the ADE reduces to a hyperbolic PDE, similar to the wave
% equation. Hyperbolic equations are characterised 
% by having both time and space derivatives of the same order.
% Hyperbolic equations can be solved analytically using the so-called
% "Method of Characteristics" (MOC).
%
% When both $a$ and $b$ are non-zero, the ADE cannot be solved using MOC and 
% it is also challenging to solve using finite differences. Situations when 
% both $a$ and $b$ are non-zero frequently come up in Earth Science. This is a
% particular problem when seeking to solve diffusion type equations, which
% have spatially varying diffusion coefficients.  For example, consider the situation:
%
% $$ \displaystyle \frac{\partial c}{\partial t}=\frac{\partial}{\partial x}
% \left(E \frac{\partial c}{\partial x}\right)=\frac{\partial E}{\partial
% x}\frac{\partial c}{\partial x}+E\frac{\partial^2 c}{\partial x^2} $$
%
% In this session we
% will examine why such difficulties occur and explore how some of the
% consequences can be dealt with.
%
% For the next few examples, the following initial and boundary
% conditions are considered:
% 
% $$ \begin{array}{lll}
% u=u_I, & 0< x \leq L, & t=0\\
% u=u_0, & x=0, & t>0\\
% u=u_I, & x=L, & t>0\\
% \end{array} $$
%  
% where $u_I$ and $u_0$ represent  specified values of $u$ and $L$ is
% the length of the domain.
%
% For the special case where $L\rightarrow\infty$, the above
% system of equations can be solved exactly to obtain
% 
% $$ \displaystyle \frac{u-u_I}{u_0-u_I}=\frac{1}{2}\left[
% \mathrm{erfc}\left(\frac{x-at}{2\sqrt{bt}}\right)
% +\exp\left(\frac{ax}{b}\right)
% \mathrm{erfc}\left(\frac{x+at}{2\sqrt{bt}}\right)
% \right]
% $$
%
% To evaluate the analytical solution, type in the following code in a 
% single script-file and save it in an appropriate place as "ADEStudy.m".
%
%   function ADEStudy
%   %A sequence of different subfunctions to explore stability and
%   %numerical diffusion associated with finite difference solutions of 
%   %the advection diffusion equation.
%
%   %Define model parameters
%   L=2; %m
%   a=1; %m/day
%   b=0.01; %m^2/day
%   uI=0; %mg/l
%   u0=1; %mg/l
%  
%   %Define location of points for plotting
%   x=linspace(0,L,50);
%   %Define values of time to be used for plotting
%   t=[1];
%   
%   %Evaluate the exact solution
%   uA=AnalSol(x,t,a,b,uI,u0);
%
%   %Plot results
%   figure(1)
%   hold off
%   plot(x,uA,'ko')
%   xlabel('Distance (m)')
%   ylabel('Concentration (mg/l)')
%   title('Concentration profile after one day')
%   
%   %**********************************************************************
%    
%   function u=AnalSol(x,t,a,b,uI,u0)
%   %Subfunction containing the analytical solution
%   [x,t]=ndgrid(x,t);
%   Pe=a*x/b;
%   F1=erfc((x-a*t)./2./sqrt(b*t))/2;
%   F2=exp(Pe).*erfc((x+a*t)./2./sqrt(b*t))/2;
%   F=F1+F2;
%   ind=Pe>700;
%   F(ind)=F1(ind);
%   u=(u0-uI)*F+uI;
%
% The outline of the remainder of this session is as follows. First we
% will derive a set of finite difference approximations using Taylor
% expansions. Following on from this we will develop four different finite
% difference solutions to the ADE along with associated stability
% criteria. We will then use the Taylor expansions to derive expressions for the
% numerical diffusion associated with each of the finite difference
% schemes. Finally, we will develop a higher order time integration scheme
% using |ode45| and compare this with the earlier first-order explicit and
% implicit time-stepping schemes.
%
%% Taylor expansions and finite difference
%
% Consider a set of discrete points in time and space, $(x_i,t_n)$, 
% where $x_i$ and $t_n$ represent the i-th point in space and n-th point 
% in time, respectively. It is assumed that these discrete points are separated by
% uniform intervals in both space and time of $\Delta x$
% and $\Delta t$, respectively. The values of a function, $u(x,t)$, at the
% locations of these discrete points can be written as $u_{i,n}= u(x_i,t_n)$.
%
% It is possible to derive finite difference approximations for the
% derivatives of $u$ with respect to $x$ and $t$ by considering the Taylor
% expansion:
%
% $$ \displaystyle u(t+\Delta t) = u(t) + \Delta t \frac{\partial u}{\partial t}+\frac{\Delta
% t^2}{2!}\frac{\partial ^2 u}{\partial t^2}+\frac{\Delta t^3}{3!}\frac{\partial ^3
% u}{\partial t^3}+O(\Delta t^4) $$
%
% from which it can be seen that:
%
% $$ \displaystyle u_{i,n+1} = u_{i,n} + \Delta t \left.\frac{\partial u}{\partial t}\right|_{i,n}
% +\frac{\Delta t^2}{2!}\left.\frac{\partial ^2 u}{\partial t^2}\right|_{i,n}
% +\frac{\Delta t^3}{3!}\left.\frac{\partial ^3
% u}{\partial t^3}\right|_{i,n}+O(\Delta t^4) $$
% 
% $$ \displaystyle u_{i+1,n} = u_{i,n} + \Delta x \left.\frac{\partial u}{\partial x}\right|_{i,n}
% +\frac{\Delta x^2}{2!}\left.\frac{\partial ^2 u}{\partial x^2}\right|_{i,n}
% +\frac{\Delta x^3}{3!}\left.\frac{\partial ^3 u}{\partial x^3}\right|_{i,n}
% +O(\Delta x^4) $$
%
% It is also possible to state that
% 
% $$ \displaystyle u(t-\Delta t) = u(t) - \Delta t \frac{\partial u}{\partial t}+\frac{\Delta
% t^2}{2!}\frac{\partial ^2 u}{\partial t^2}-\frac{\Delta t^3}{3!}\frac{\partial ^3
% u}{\partial t^3}+O(\Delta t^4) $$
% 
% from which it can be seen that:
% 
% $$ \displaystyle u_{i,n} = u_{i,n+1} - \Delta t \left.\frac{\partial u}{\partial t}\right|_{i,n+1}+\frac{\Delta
% t^2}{2!}\left.\frac{\partial ^2 u}{\partial t^2}\right|_{i,n+1}-\frac{\Delta t^3}{3!}\left.\frac{\partial ^3
% u}{\partial t^3}\right|_{i,n+1}+O(\Delta t^4) $$
% 
% $$ \displaystyle u_{i-1,n} = u_{i,n} - \Delta x \left.\frac{\partial u}{\partial x}\right|_{i,n}+\frac{\Delta
% x^2}{2!}\left.\frac{\partial ^2 u}{\partial x^2}\right|_{i,n}-\frac{\Delta x^3}{3!}\left.\frac{\partial ^3
% u}{\partial x^3}\right|_{i,n}+O(\Delta x^4) $$
% 
% Interestingly, it can also be seen that:
% 
% $$ \displaystyle u_{i+1,n}-u_{i-1,n} =2 \Delta x \left.\frac{\partial u}{\partial x}\right|_{i,n}+2\frac{\Delta x^3}{3!}\left.\frac{\partial ^3
% u}{\partial x^3}\right|_{i,n}+O(\Delta x^4) $$
% 
% $$ \displaystyle u_{i+1,n}+u_{i-1,n} = 2u_{i,n}+2\frac{\Delta
% x^2}{2!}\left.\frac{\partial ^2 u}{\partial x^2}\right|_{i,n}+O(\Delta x^4)
% $$
% 
% Rearranging the above set of equations leads to
% the following set of finite difference approximations:
% 
% $$ \displaystyle \frac{u_{i+1,n}-2u_{i,n}+u_{i-1,n}}{\Delta x^2}
% =\left.\frac{\partial ^2 u}{\partial x^2}\right|_{i,n}+O(\Delta
% x^2) 
% $$ _____(2a)
% 
% $$ \displaystyle \frac{u_{i+1,n}-u_{i-1,n}}{2 \Delta x}= \left.\frac{\partial
% u}{\partial x}\right|_{i,n} +O(\Delta x^2) 
% $$ _____(2b)
% 
% $$ \displaystyle \frac{u_{i,n}-u_{i-1,n}}{\Delta x} = \left.\frac{\partial
% u}{\partial x}\right|_{i,n} -\frac{\Delta
% x}{2!}\left.\frac{\partial ^2 u}{\partial x^2}\right|_{i,n}
% +O(\Delta x^2) $$ _____(2c)
% 
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t} =  \left.\frac{\partial
% u}{\partial t}\right|_{i,n+1} -\frac{\Delta
% t}{2!}\left.\frac{\partial ^2 u}{\partial t^2}\right|_{i,n+1}
% +O(\Delta t^2) $$ _____(2d)
% 
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t} =  \left.\frac{\partial
% u}{\partial t}\right|_{i,n} +\frac{\Delta
% t}{2!}\left.\frac{\partial ^2 u}{\partial t^2}\right|_{i,n}
% +O(\Delta t^2) $$ _____(2e)
%
%% ETCS - Explicit time-stepping with central differences in space
%
% First we will develop an Euler explicit time-stepping scheme for the ADE 
% problem described above and use a
% central difference approximation for the $\partial u/\partial x$ term.
% Ignoring all the truncation errors associated with the Taylor expansions,
% the ADE can then be written as
%
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t} +a\left(\frac{u_{i+1,n}-u_{i-1,n}}
% {2 \Delta x}\right)
% -b\left(\frac{u_{i+1,n}-2u_{i,n}+u_{i-1,n}}{\Delta x^2}\right)=0 $$ 
% _____(3a)
%
% which can be rearranged to get
%
% $$ u_{i,n+1}=A u_{i-1,n}+B u_{i,n} + C u_{i+1,n} $$
% 
% where
% 
% $$ \displaystyle A=\frac{a\Delta t}{2 \Delta x}+\frac{b\Delta t}{\Delta x^2},
% \quad
% B=1-\frac{2b\Delta t}{\Delta x^2},
% \quad
% C=-\frac{a\Delta t}{2 \Delta x}+\frac{b\Delta t}{\Delta x^2} $$
% 
% Hereafter this scheme will be referred to as the ETCS scheme.
%
% For numerical stability of explicit time-stepping schemes to be ensured it is
% generally found that the coefficients, $A$, $B$ and $C$, need
% always to be positive. Given that $a$, $b$, $\Delta x$ and $\Delta
% t$ should also always be positive, it follows that there are two
% stability criteria for this scheme:
% 
% $$ \displaystyle \mathrm{Cr}\equiv \frac{b\Delta t}{\Delta x^2}<0.5 
% \quad
% \mathrm{and}
% \quad
% \mathrm{Pe}\equiv  \frac{a\Delta x}{b }<2 $$
% 
% where $\mathrm{Cr}$ and $\mathrm{Pe}$ are often referred
% to as the Courant number and Peclet number, respectively.
%
% To evaluate the ETCS scheme add the following subfunction to ADEStudy.m:
%
%   %**********************************************************************
%
%   function u=FDSol(xI,tI,a,b,uI,u0,Pe,Cr)
%   %Subfunction containing the finite difference solutions
%
%   %Determine the space-step and time-step using a specified
%   %Courant and Peclet number
%   dx=Pe*b/a;
%   dt=Cr*dx^2/b;
%   
%   %Determine number of nodes to solve for
%   Nx=ceil((xI(end)-xI(1))/dx)+1;
%   Nt=ceil(tI(end)/dt)+1;
%   
%   %Determine the location of solution points in x and t
%   x=[0:dx:dx*(Nx-1)]+xI(1);
%   t=[0:dt:dt*(Nt-1)];
%   
%   %Initialise the solution matrix
%   u=zeros(Nx,Nt);
%   
%   %Apply initial condition
%   u(:,1)=uI;
%   
%   %Apply boundary conditions
%   u(1,:)=u0;
%   u(Nx,:)=uI;
%   
%   %Determine A, B and C coefficients
%   A=a*dt/2/dx+b*dt/dx^2;
%   B=1-2*b*dt/dx^2;
%   C=-a*dt/2/dx+b*dt/dx^2;
%   
%   %Determine Jacobian matrix
%   J=spdiags(ones(Nx,1)*[A B C],[-1 0 1],Nx,Nx);
%
%   %Apply boundary conditions
%   J(1,:)=0;
%   J(1,1)=1;
%   J(Nx,:)=0;
%   J(Nx,Nx)=1;
%   
%   %Solve problem
%   for n=1:Nt-1
%       u(:,n+1)=J*u(:,n);
%   end
%   
%   %Interpolate solution to plotting points
%   u=interp2(t,x',u,tI,xI');
%
% Note that the values of $\Delta x$ and $\Delta t$ are determined by
% considering a specified Peclet and Courant number. Also note that the
% solution is then interpolated to the values of $x$ and $t$ specified by
% the input arguments of the sub-function. This removes the need to plot
% all the solution points solved for, which can be excessive depending on
% the Peclet number and Courant number specified.
%
% It is also interesting to see that the $A$, $B$ and $C$ coefficients represent
% the values of the first, second and third diagonal of the 
% Jacobian matrix, $\mathbf{J}$, defined in this case by the formula
%
% $$ \mathbf{u}_{n+1}=\mathbf{J}\mathbf{u}_{n} $$
%
% To evaluate the ETCS scheme, replace:
%
%   %Plot results
%   figure(1)
%   hold off
%   plot(x,uA,'k')
%   xlabel('Distance (m)')
%   ylabel('Concentration (mg/l)')
%   title('Concentration profile after one day')
%
% in the main function of ADEStudy.m with:
%
%   %Specify a Peclet number and Courant number
%   Pe=2;
%   Cr=0.5;
%   
%   %Evaluate the finite difference solution
%   uFD=FDSol(x,t,a,b,uI,u0,Pe,Cr);
%   
%   %Plot results
%   figure(1)
%   hold off
%   plot(x,uA,'k')
%   hold on
%   plot(x,uFD)
%   xlabel('Distance (m)')
%   ylabel('Concentration (mg/l)')
%   title('Concentration profile after one day')
%   legend('Analytical solution','Finite difference solution')
%
% Run the script-file and compare the results from the analytical solution
% and the finite difference solution. Experiment by varying the specified
% value of $\mathrm{Pe}$ and $\mathrm{Cr}$. Do the pre-determined stability 
% criteria work? Does stability also ensure accuracy?
%
%% ITCS - Implicit time-stepping with central differences in space
%
% Here we will develop an Euler implicit time-stepping scheme and use a
% central difference approximation for the $\partial u/\partial x$ term.
% Ignoring all the truncation errors associated with the Taylor expansions,
% the ADE can then be written as
%
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t}
% +a\left(\frac{u_{i+1,n+1}-u_{i-1,n+1}}{2 \Delta
% x}\right)-b\left(\frac{u_{i+1,n+1}-2u_{i,n+1}+u_{i-1,n+1}}{\Delta
% x^2}\right)=0 $$ _____(3b)
%
% which can be rearranged to get
%
% $$ u_{i,n}=A u_{i-1,n+1}+B u_{i,n+1} + C u_{i+1,n+1} $$
% 
% where
% 
% $$ \displaystyle A=-\frac{a\Delta t}{2 \Delta x}-\frac{b\Delta t}{\Delta x^2},
% \quad
% B=1+\frac{2b\Delta t}{\Delta x^2},
% \quad
% C=\frac{a\Delta t}{2 \Delta x}-\frac{b\Delta t}{\Delta x^2} $$
% 
% Hereafter this scheme will be referred to as the ITCS scheme.
%
% For numerical stability of implicit time-stepping schemes to
% be ensured it is
% generally found that the coefficients, $A$ and $C$ need
% always to be negative whereas $B$ needs to be positive.
% Given that $a$, $b$, $\Delta x$ and $\Delta
% t$ should always be positive, it follows that there  remains still one
% stability criterion:
% 
% $$ \displaystyle \mathrm{Pe}\equiv  \frac{a\Delta x}{b }<2 $$
% 
%
% To evaluate and compare the ITCS scheme to the ETCS scheme and the
% analytical solution we need to make the following modifications to ADEStudy.m.
%
% First we need to add a term called |Scheme| as an input argument to the 
% |FDSol| subfunction. This can be done by modifying the start of the
% subfunction to read as follows.
%
%   function u=FDSol(xI,tI,a,b,uI,u0,Pe,Cr,Scheme)
%
% The term, |Scheme|, can be set to |'ITCS'| or |'ETCS'|. This in turn
% can be used to dictate which expressions should be used to calculate 
% the $A$, $B$ and
% $C$ coefficients. To enable this flexibility replace the code:
%
%   %Determine A, B and C coefficients
%   A=a*dt/2/dx+b*dt/dx^2;
%   B=1-2*b*dt/dx^2;
%   C=-a*dt/2/dx+b*dt/dx^2;
%
% in the  |FDSol| subfunction with:
%
%   %Determine A, B and C coefficients
%   switch Scheme
%       case 'ETCS'
%           %Explicit time-stepping central difference in space
%           A=a*dt/2/dx+b*dt/dx^2;
%           B=1-2*b*dt/dx^2;
%           C=-a*dt/2/dx+b*dt/dx^2;
%       case 'ITCS'
%           %Implicit time-stepping central difference in space
%           A=-a*dt/2/dx-b*dt/dx^2;
%           B=1+2*b*dt/dx^2;
%           C=a*dt/2/dx-b*dt/dx^2;
%   end
%
% The above modification allows the use of |Scheme| to
% determine whether the user requires the ETCS or the ITCS scheme. Based on
% that information, the |switch| function can select the appropiate set of
% expressions for the $A$, $B$ and $C$ coefficients. Read about |switch| in the help
% documentation to learn more about |switch| and how it can be used.
%
% Again it should be noted that the $A$, $B$ and $C$ coefficients represent
% the values of the first, second and third diagonals of a tridiagonal
% Jacobian matrix, $\mathbf{J}$. But for the ITCS scheme, $\mathbf{J}$ is 
% defined by the formula
%
% $$ \mathbf{u}_{n}=\mathbf{J}\mathbf{u}_{n+1} $$
%
% MATLAB's left divide operator can be used  to solve such a system of linear
% equations. For example, |x=A\B| solves the system of linear equations
% |A*x=B|.
%
% To implement the above equation replace
%
%   %Solve problem
%   for n=1:Nt-1
%       u(:,n+1)=J*u(:,n);
%   end
%
% with
%
%   %Solve problem
%   switch Scheme(1:2)
%       case 'ET'
%           %Explicit time-stepping
%           for n=1:Nt-1
%               u(:,n+1)=J*u(:,n);
%           end
%       case 'IT'
%           %Implicit time-stepping
%           for n=1:Nt-1
%               u(:,n+1)=J\u(:,n);
%           end
%   end
%
% Note that we only need to use the first two characters of the string contained
% in |Scheme|, because the above code is only concerned with whether the scheme
% involves implicit time-stepping (i.e., |'IT'|) or explicit time-stepping
% (i.e., |'ET'|).
%
% To compare the analytical solution with ETCS and ITCS schemes
% collectively, replace the following code within the main function
%
%   %Evaluate the finite difference solution
%   uFD=FDSol(x,t,a,b,uI,u0,Pe,Cr);
%
% with
%
%   %Evaluate the explicit time-stepping with central differences in space scheme
%   uFD(:,1)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ETCS');
%   %Evaluate the implicit time-stepping with central differences in space scheme
%   uFD(:,2)=FDSol(x,t,a,b,uI,u0,Pe,Cr,'ITCS');
%
% and update the legend appropriately.
%
% Run the modified script-file and compare the results from the analytical solution
% and the finite difference solutions. Experiment by varying the specified
% value of $\mathrm{Pe}$ and $\mathrm{Cr}$.
%
%% ITBS - Implicit time-stepping with backward differences
%
% If we consider an Euler implicit time-stepping scheme and use a
% backward difference approximation for the $\partial u/\partial x$ term,
% ignoring all the truncation errors associated with the Taylor expansions,
% the ADE can then be written as
%
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t}
% +a\left(\frac{u_{i,n+1}-u_{i-1,n+1}}{ \Delta
% x}\right)-b\left(\frac{u_{i+1,n+1}-2u_{i,n+1}+u_{i-1,n+1}}{\Delta
% x^2}\right)=0 $$ _____(3c)
%
% which can be rearranged to get
%
% $$ u_{i,n}=A u_{i-1,n+1}+B u_{i,n+1} + C u_{i+1,n+1} $$
% 
% where
% 
% $$ \displaystyle A=-\frac{a\Delta t}{ \Delta x}-\frac{b\Delta t}{\Delta x^2},
% \quad
% B=1+\frac{a\Delta t}{ \Delta x}+\frac{2b\Delta t}{\Delta x^2},
% \quad
% C=-\frac{b\Delta t}{\Delta x^2} $$
%
% from which it can be seen that this scheme is unconditionally stable.
%
% Hereafter this scheme will be referred to as the ITBS scheme.
%
%% ETBS - Explicit time-stepping with backward differences
%
% If we consider an Euler explicit time-stepping scheme and use a
% backward difference approximation for the $\partial u/\partial x$ term,
% ignoring all the truncation errors associated with the Taylor expansions,
% the ADE can then be written as
%
% $$ \displaystyle \frac{u_{i,n+1}-u_{i,n}}{ \Delta t}
% +a\left(\frac{u_{i,n}-u_{i-1,n}}{ \Delta
% x}\right)-b\left(\frac{u_{i+1,n}-2u_{i,n}+u_{i-1,n}}{\Delta
% x^2}\right)=0 $$ _____(3d)
%
% which can be rearranged to get
%
% $$ u_{i,n+1}=A u_{i-1,n}+B u_{i,n} + C u_{i+1,n} $$
% 
% where
% 
% $$ \displaystyle A=\frac{a\Delta t}{ \Delta x}+\frac{b\Delta t}{\Delta x^2},
% \quad
% B=1-\frac{a\Delta t}{ \Delta x}-\frac{2b\Delta t}{\Delta x^2},
% \quad
% C=\frac{b\Delta t}{\Delta x^2} $$
%
% from which it can be seen that this scheme is conditionally stable
% providing
%
% $$ \displaystyle  \mathrm{Cr}<\frac{1}{2+\mathrm{Pe}} $$
%
% Hereafter this scheme will be referred to as the ETBS scheme.
%
%% Classroom Assignment - Part 1
%
% Modify the ADEStudy.m file further to also include and compare the ITBS
% and ETBS schemes. Generate a single figure overlaying the plots of
% concentration against distance using each of the different schemes. How
% do they compare with the results from the analytical solution?
%
%% Numerical diffusion
%
% It is apparent that the stability criteria for the different schemes
% listed above do not necessarily to lead to accurate solutions. If you set
% $\mathrm{Pe}=2$ and $\mathrm{Cr}=1/(2+\mathrm{Pe})$, you will see that
% all the finite difference solutions overestimate the amount of diffusion
% taking place, with the exception of the ETCS scheme.
%
% It is possible to determine exactly how much the different schemes are
% overestimating the diffusion by consideration of the truncation errors in
% the Taylor expansions provided earlier.
%
% _The ETCS scheme_
%
% Consider again the ETCS scheme (i.e., Eq. (3a)). Substituting Eqs. (2a),
% (2b) and (2e) into Eq. (3a) whilst retaining the truncation errors leads
% to
%
% $$ \displaystyle \frac{\partial u}{\partial t}
% +\frac{\Delta t}{2}\frac{\partial ^2 u}{\partial t^2}
% + a\frac{\partial u}{\partial x}-
% b\frac{\partial ^2 u}{\partial x^2}
% =O(\Delta x^2,\Delta t^2) $$ _____(4a)
%
% Now consider Eq. (1), it follows that
%
% $$ \displaystyle \frac{\partial ^2 u}{\partial t^2}=\frac{\partial}{\partial t}
% \left(-a\frac{\partial u}{\partial x}+b\frac{\partial ^2 u}{\partial x^2}\right) $$
% 
% After some further rearranging (and recalling that $a$ and $b$ are constants)
% it can be further shown that
% 
% $$ \displaystyle \frac{\partial ^2 u}{\partial t^2}=a^2\frac{\partial ^2 u}{\partial x^2}
% -2ab\frac{\partial ^3 u}{\partial x^3}+b^2 \frac{\partial ^4 u}{\partial
% x^4} $$ _____(4b)
% 
% which on substitution into Eq. (4a) leads to
% 
% $$ \displaystyle \frac{\partial u}{\partial t}
% + a\frac{\partial u}{\partial x}-
% \left(b-\frac{a^2\Delta t }{2}\right)\frac{\partial ^2 u}{\partial x^2}
% =ab\Delta t\frac{\partial ^3 u}{\partial x^3}
% -\frac{b^2 \Delta t}{2}\frac{\partial ^4 u}{\partial x^4}+O(\Delta x^2,\Delta t^2) $$
% 
% from which it can be seen that the ETCS scheme gives rise to a
% numerical diffusion, $b_{num}$, of magnitude
%
% $$ \displaystyle b_{num}=-\frac{a^2\Delta t }{2}=-\frac{b \mathrm{Pe}^2\mathrm{Cr}}{2}
% $$
%
% which is negative, hence explaining why ETCS appears to
% underestimate the amount of diffusion predicted by the analytical
% solution.
%
% _The ITCS scheme_
%
% Consider again the ITCS scheme (i.e., Eq. (3b)). Substituting Eqs. (2a),
% (2b) and (2d) into Eq. (3b) and then substituting Eq. (4b), whilst 
% retaining the truncation errors, leads to
%
% $$ \displaystyle \frac{\partial u}{\partial t}
% + a\frac{\partial u}{\partial x}-
% \left(b+\frac{a^2\Delta t }{2}\right)\frac{\partial ^2 u}{\partial x^2}
% =\frac{b^2 \Delta t}{2}\frac{\partial ^4 u}{\partial x^4}-ab\Delta t\frac{\partial ^3 u}{\partial x^3}
% +O(\Delta x^2,\Delta t^2) $$
%
% from which it can be seen that the ITCS scheme gives rise to a
% numerical diffusion, $b_{num}$, of magnitude
%
% $$ \displaystyle b_{num}=\frac{a^2\Delta t }{2}=\frac{b \mathrm{Pe}^2\mathrm{Cr}}{2}
% $$
%
% _The ITBS scheme_
%
% Consider again the ITBS scheme (i.e., Eq. (3c)). Substituting Eqs. (2a),
% (2c) and (2d) into Eq. (3c) and then substituting Eq. (4b), whilst 
% retaining the truncation errors, leads to
%
% $$ \displaystyle \frac{\partial u}{\partial t}
% + a\frac{\partial u}{\partial x}-
% \left(b+\frac{a\Delta x}{2}+\frac{a^2\Delta t }{2}\right)\frac{\partial ^2 u}{\partial x^2}
% =\frac{b^2 \Delta t}{2}\frac{\partial ^4 u}{\partial x^4}-ab\Delta t\frac{\partial ^3 u}{\partial x^3}
% +O(\Delta x^2,\Delta t^2) $$
%
% from which it can be seen that the ITBS scheme gives rise to a
% numerical diffusion, $b_{num}$, of magnitude
%
% $$ \displaystyle b_{num}=\frac{a\Delta x}{2}+\frac{a^2\Delta t }{2}=\frac{b \mathrm{Pe}(1+\mathrm{Pe}\mathrm{Cr})}{2}
% $$
%
% _The ETBS scheme_
%
% Consider again the ETBS scheme (i.e., Eq. (3d)). Substituting Eqs. (2a),
% (2c) and (2e) into Eq. (3d) and then substituting Eq. (4b), whilst 
% retaining the truncation errors, leads to
%
% $$ \displaystyle \frac{\partial u}{\partial t}
% + a\frac{\partial u}{\partial x}-
% \left(b+\frac{a\Delta x}{2}-\frac{a^2\Delta t }{2}\right)\frac{\partial ^2 u}{\partial x^2}
% =ab\Delta t\frac{\partial ^3 u}{\partial x^3}
% -\frac{b^2 \Delta t}{2}\frac{\partial ^4 u}{\partial x^4}+O(\Delta x^2,\Delta t^2) $$
%
% from which it can be seen that the ETBS scheme gives rise to a
% numerical diffusion, $b_{num}$, of magnitude
%
% $$ \displaystyle b_{num}=\frac{a\Delta x}{2}-\frac{a^2\Delta t }{2}=\frac{b \mathrm{Pe}(1-\mathrm{Pe}\mathrm{Cr})}{2}
% $$
%
%% Classroom Assignment - Part 2
%
% 1) Modify the ADEStudy.m file further to also evaluate the analytical
% solution with $b=b+b_{num}$ using the four different values of $b_{num}$
% associated with the four finite difference schemes discussed above. Plot
% the results as dot markers over the finite difference results. How
% effective are the expressions for $b_{num}$ at predicting the diffusion
% associated with each of the finite difference schemes? Which finite
% difference scheme do you think is best and why?
%
% 2) Solve the above problem using |ode45| for the time integration and try
% using central differences and then backward differences for the $\partial
% u/\partial x$ term. Compare your results with those from the other finite
% difference schemes. Also, derive expressions for the numerical diffusion
% associated with the two new schemes and check these work by evaluating
% the analytical solution with the additional diffusion terms (as in part 2
% of the classroom assignment above).
