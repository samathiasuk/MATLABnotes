%% MATLAB Notes: Session 4 - Approximate methods for differentiation and integration
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
% * Describe how derivatives can be approximated using finite differences.
% * Approximate the derivative of a variable using |diff|.
% * Show how finite differences link with trapezoidal integration.
% * Approximate the integral of a variable using |cumtrapz|.
%
%% Finite differences
%
%
% The basic idea about finite differences is as follows:
%
% Consider a continuous function $y=f(x)$.
%
% Now consider a set of $N$ number of discrete points along the $x$ 
% axis $x=x_1,x_2,x_3,...,x_N$.
%
% The corresponding values of $y$ can be written as $y_1,y_2,y_3,...,y_N$.
%
% Understand that $y(x=x_1)=y_1$, $y(x=x_2)=y_2$ and so on.
%
% An approximation of the derivative
% of $y$ with respect to $x$ can be obtained as follows:
%
% $\displaystyle \left. \frac{dy}{dx} \right|_{x=x_{i-1/2}} \approx
% \frac{y_i-y_{i-1}}{x_i-x_{i-1}}$ where $i = 2,3,...,N$.     
%
% The value of the derivative corresponds to the location, $x=x_{i-1/2}$
% where $x_{i-1/2}=(x_i+x_{i-1})/2$.
%
% Note that calculus is often referred to as "infinitesimal
% calculus". The infinitesimal term refers to the focus on how $y$ 
% changes over infinitesimal distances. A definition of a
% derivative can be written as follows:
%
% $\displaystyle \lim_{x_i\rightarrow x_{i-1}}\left. \frac{dy}{dx}
% \right|_{x=x_{i-1/2}} =
% \frac{y_i-y_{i-1}}{x_i-x_{i-1}}$ 
% 
%% Approximating derivatives using finite differences
%
% Consider the trigonometric equation
%
% $y = \cos 2x$
%
% Differentiating with respect to $x$ gives
%
% $\displaystyle \frac{dy}{dx} = -2\sin 2x$
%
% Creat a new script file and type:
 
% Some examples demonstrating the principal of differentiation and 
% integration using finite  differences and trapezoidal integration, 
% respectively.
 
%%
% at the top and save in an appropriate place.
%
% First we will make a series of $N$ number of equally spaced points from zero to five
% using the |linspace| command (type "help linspace" in the Command Window
% to find out more). Following from this will calculate some corresponding
% values of $y$ (where $y=\cos 2x$).
 
%Define number of points for the discretisation of x
N=7;
x=linspace(0,5,N)
%Calculate corresponding values of y
y=cos(2*x)
 
%%
% Now we will estimate the derivative using finite differences. Note that
% you can use |end| in the code to automatically identify the number of
% elements for a given dimension of an array of concern.
 
%Calculate the change in x, i.e. x(i)-x(i-1)
dx=x(1,2:end)-x(1,1:end-1)
%Calculate the change in y, i.e. y(i)-y(i-1)
dy=y(1,2:end)-y(1,1:end-1)
%Estimate the derivative with respect to x using finite difference
dydxFD=dy./dx
 
%%
% Interestingly, we have only six values of the derivative although we started with seven
% points. It is important to realise that these derivatives do not
% correspond directly to the points:
%
% $x_1,x_2,x_3,...,x_N$
%
% Rather they correspond to a different set of points:
%
% $x_{2-1/2},x_{3-1/2},x_{4-1/2},...,x_{N-1/2}$
%
% of which there are only $(N-1)$ in number.
%
% As discussed earlier, it is possible to locate this latter set of points from
%
% $x_{i-1/2}=(x_i+x_{i-1})/2$, where $i=2,3,...,N-1$
%
% Add the following to your code:
 
%Calculate the locations at which the derivatives apply
xMP=(x(1,2:end)+x(1,1:end-1))/2
 
%%
% Recall that we know the derivative is, in this case, 
% $-2\sin 2x$. So lets compare:
 
%Calculate derivative analytically at the xMP points
dydz=-2*sin(2*xMP)
 
%%
% Now lets make a plot to compare our results graphically
 
%Plot results
figure(1)
clf
plot(x,y,'k',xMP,dydz,'b',xMP,dydxFD,'ro')
legend('y = cos 2x','dy/dx Analytical','dy/dx Finite difference')
xlabel('x (-)')
ylabel('see legend')
title('Approximate calculus')
 
%%
% The accuracy is quite poor because we are using only seven
% points for our discretisation. Try running your code again with 20
% points.
%
% Although using finite difference techniques to differentiate trigonometric
% functions seems slightly trivial, we are just using these as a
% demonstration. The main point to understand is that we can use similar
% code to approximate the derivative of any function of interest.
%
%% Approximating derivatives using the diff command
%
% The |diff| command is a MATLAB command, which can be used to perform the
% above analysis in a computationally more efficient manner. Type "help
% diff" in the command window to find out more.
%
% Note that |diff(X,N,DIM)|, where |X| is the name of the array to be differenced, |N| is the number
% of times the array us to be differences, and |DIM| is the dimension of
% the array along which the differencing is to take place.
%
% To aid our demonstration of |diff|, type the following 2D array in the command window:
 
A=[1 6 5 4 3 2; 5 6 4 4 2 1; 7 8 5 1 1 2]
 
%%
% Now we are going to look at the difference between each column
 
A(:,2:end)-A(:,1:end-1)
 
%%
% The |diff| command can be used to do exactly the same thing by typing:
 
diff(A,1,2)
 
%%
% If instead we wanted to look at the difference between each row, we would
% type:
 
diff(A,1,1)
 
%%
% If we wanted to look at the difference of each difference between each
% row, we would type:
 
diff(A,2,1)
 
%%
% Go back to your script file that you created and add the following code
 
dydxFD=diff(y,1,2)./diff(x,1,2)
 
%%
% Above we have calculated the derivative in one line. Note that the |1|, 
% in the |diff| argument, is
% because we are only taking one difference and the |2| is because we want
% to difference the columns as opposed to the rows. Recall that |x| and
% |y| are row vectors.
 
%% Trapezoidal rule and its relation to finite difference
%
% Recall that integration is the reverse of differentiation. Whereas
% differentiation corresponds to finding the gradient of a line on a graph,
% integration corresponds to finding the area under a line.
%
% A simple way to approximate an integral is to apply the so-called
% trapezoidal rule. Click on this  link
% http://mathworld.wolfram.com/TrapezoidalRule.html
% TrapezoidalRule for an illustration.
%
% Recalling that integration is the reverse of differentiation. Consider
% again the finite difference approximation
%
% $\displaystyle \left. \frac{dy}{dx} \right|_{x=x_{i-1/2}} \approx
% \frac{y_i-y_{i-1}}{x_i-x_{i-1}}$ where $i = 2,3,...,N$
%
% For convenience, let us say $dy/dx=J$ such that
%
% $\displaystyle J_{i-1/2} \approx
% \frac{y_i-y_{i-1}}{x_i-x_{i-1}}$
%
% It can also be said
%
% $\displaystyle y=\int J dx$
%
% Rearranging the finite difference equation above such that $y_i$ is the
% subject of the formula we get
%
% $y_i \approx J_{i-1/2}(x_i-x_{i-1})+y_{i-1}$
%
% In the same way that $x_{i-1/2}=(x_i+x_{i-1})/2$, it can also be said
% that $J_{i-1/2}=(J_i+J_{i-1})/2$, from which it can be seen 
%
% $\displaystyle y_i \approx \left(\frac{J_{i}+J_{i-1}}{2}\right)(x_i-x_{i-1})+y_{i-1}$
%
% which is essentially a mathematical expression of the trapezoidal rule.
 
%% Approximating integrals using the trapezoidal rule
%
% Consider again the trigonometric equation
%
% $y = \cos 2x$
%
% Integrating with respect to $x$ gives
%
% $\displaystyle F = \int y dx = \frac{1}{2}\sin 2x + C$
% 
% where $C$ is an integration constant.
%
% Let us now impose the constraint that $F=3$ when $x =0$. It
% follows that $C = 3$. Therefore 
%
% $\displaystyle F = 3+\frac{1}{2}\sin 2 x$
%
% We now wish to approximate $F$ using the trapezoidal rule.
%
% From the previous discussion we can say that
%
% $\displaystyle F_i \approx \left(\frac{y_{i}+y_{i-1}}{2}\right)(x_i-x_{i-1})+F_{i-1}$
%
% To explore this further let us add the following code to our script file:
 
%Initialise a vector of zeros called Ftrapz
Ftrapz=zeros(1,N);
%Estimate the integral of y with respect to x using the trapezoidal rule
for i=2:N
    Ftrapz(1,i)=(y(1,i)+y(1,i-1))/2*(x(1,i)-x(1,i-1))+Ftrapz(1,i-1);
end
%Display result in the command window
disp(Ftrapz)
 
%%
% Now lets add some code to compare the result to the analytical solution:
 
%Calculate values F using the analytical function from notes
F=3+sin(2*x)/2
 
%%
% Note that we also need to apply our constraint $ F(x=0)=3 $ to the trapezoidal
% integration result. This can be done as follows:
 
%Apply constraint to trapezoidal integration
Ftrapz=Ftrapz-Ftrapz(1,1)+3
 
%%
% What we have done is subtracted the first value of |Ftrapz|, which 
% corresponds to $x = 0$ in this case, from all the values of |F|. We
% have then added on the correct value at $x = 0$, which according to
% our constraint, is 3.
 
%Plot results
figure(1)
hold on %Plot new results in figure(1) alongside previous results
plot(x,F,'g',x,Ftrapz,'k:s')
legend('y = cos 2x','dy/dx Analytical','dy/dx Finite difference',...
    '\int y dx Analytical','\int y dx Trapezoidal')
legend('location','eastoutside')
 
%% Approximating integrals using the cumtrapz command
%
% The |cumtrapz| command is a MATLAB command, which can be used to perform
% exactly the same analysis above in a computationally more efficient manner. Type "help
% cumtrapz" in the command window to find out more.
%
% Note that |cumtrapz(X,Y,DIM)|, where |X| is the name of the array to be 
% integrated with respect to, |Y| is the array to be integrated, and |DIM| is the dimension of
% the array along which the integrating is to take place.
%
% Add the following code to your script file to repeat the integration
% using |cumtrapz|:
 
%Use cumtrapz to integrate y with respect to x
Fcumtrapz=cumtrapz(x,y,2);
%Apply constraint
Fcumtrapz=Fcumtrapz-Fcumtrapz(1,1)+3
 
%% Classroom assignment
%
% Make a new script file and write a subfunction to both differentiate and
% integrate a given set of dependent and independent variables using |diff| 
% and |cumtrapz| (if one differentiates $y$ with respect to $x$,
% $y$ and $x$ are known as the dependent and independent variables,
% respectively).
%
% Then write some additional script to use your subfunction to 
% both differentiate and integrate the following
% functions  for $0 < x < 2$:
%
% 1) $y = 5x^2 + 2x +2$
%
% 2) $y = e^x$
%
% 3) $y = (1+x^2)^{-1}$
%
% For each equation apply a constraint such that the integrals all equal 1
% when $x = 0$.
%
% Write further code to plot your results appropriately on a single
% figure with three subplots. Compare your results with analytical results
% where possible. You may find www.wolframalpha.com useful in this respect.
% 
% An example MATLAB code, containing all of the instructions requested
% above, is given in <https://github.com/samathiasuk/MATLABnotes/blob/main/MATLABsession4_Assignment.m MATLABsession4_Assignment.m>.
 
 


