%% MATLAB Notes: Session 2 - Loops and conditions
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
% * Repeat operations using a "for" loop.
% * Define logical statements.
% * Compare arrays using logical statements.
% * Select data using logical statements.
% * Use sum, mean, std, max and min.
%
%% Repeat operations using a "for" loop
%
% Last week we learnt how to store and access data in arrays. We also
% learnt how to write long lists of instructions in a script file.
% Sometimes we want to do the same or a similar thing many times. 
% Instead of repeating the code many times we can use a "for" loop.
%
% Open MATLAB.
%
% Click on "File", "New", "Script". The "Editor" window should now appear.
% Save the file as "MATLABsession2_Assignment.m" in the your "MATLABnotes"
% directory.
%
% Type the following code in your new mfile and run (by pressing the green
% arrow at the top of the editor window).

%This is a script file containing examples from MATLAB session 2
A=2;
A=A+5;
A=A+5;
A=A+5;
%Display the value of A in the Command Window
A
 
%%
% What we have done is added 5 onto 2, three times. Another way of doing
% this is to use a for loop. Add the following to your mfile and run:
 
A=2;
for i=1:3
    %Display the value of i in the Command Window
    i
    %Add 5 onto A and display in the Command Window
    A=A+5
end
 
%%
% Exactly the same result is achieved. We have added 5 onto 2, three times.
% So what exactly is happening? We are stepping through a loop. In the
% first step, |i=1|. In the second step, |i=2|, in the third step, |i=3|.
% Then the loop is stopped. During each step 5 is added onto |A|.
%
% Consider what happens when we type
 
1:3
 
%%
% At the beginning of the loop we state |for i=1:3|. This tells
% MATLAB to make |i=1|, then |2| and finally |3|. Anything written
% in between the |for| and |end| is repeated until the sequence of |i|
% values is achieved.
%
% Add the following to your mfile and run:
 
A=2;
%Add 5 on to A 1000 times
for i=1:1000
    A=A+5;
end
%Display the value of A in the Command Window
A
 
%%
% We have repeated the calculation 1000 times!
%
%% Defining logical statements
%
% Another important aspect of scientific programming concerns the application
% of logical statements.
%
% Conditions usually contain the following:
 
% ==    equal to
% <     less than
% <=    less than or equal to
% >     greater than
% >=    greater than or equal to
% ~=    not equal to
 
%%
% MATLAB denotes true as 1 and false as 0. For example, type the following
% into the Command Window:
 
5<6
 
%%
% Indeed it is "true" that 5<6, therefore MATLAB returned a 1.
%
% Type the following into the Command Window:
 
5>6
 
%%
% It is "false" that 5>6, therefore MATLAB returned a 0.
%
% Type the following into the Command Window
 
x=[1 3 4 8 6 5 3 2 2]
x<3
 
%%
% The |x<3| leads to a row of ones and zeros. There is a one for every
% value of x which is < 3 and a zero for every value of x which is greater 
% than or equal to 3.
%
% Now type: 
 
x<=3
 
%%
% 
% Now there is a one for every value of x which is less than or equal to 3 
% and a zero for every value of x which is > 3.
%
% What about if we want find out which values are equal to 3? Type

x=[1 3 4 8 6 5 3 2 2]
x=3
 
%%
% Oh dear, now x = 3. The makers of MATLAB realised this is the case and 
% so introduced the |==| term.
%
% Now type:
 
x=[1 3 4 8 6 5 3 2 2]
x==3
 
%%
% We get a one for every value of x which is equal to 3. Now type
 
x~=3
 
%%
% Now we get a one for every value of x which is not equal to 3.
%
%% Comparing arrays using logical statements.
%
% We use these statements to compare different arrays. However, 
% the arrays must be the same size. Type
 
x=[1 3 4 8 6 5 3 2 2]
y=[3 3 1 5 3 2 8 4 2]
x==y
 
%%
% We get a one for every value of x which is equal to the corresponding 
% value of y.
 
 
%%
% Another two useful symbols are
 
% &    and
% |    or
 
%%
% Suppose we want to know about which values of x are equal to y and are
% also greater than 2. Type
 
x==y & x>2
 
%%
% Suppose we want to know which values of x are greater than y or equal to
% 2. Type
 
x>y | x==2
 
%%
% We can make these statements as complicated as we like. 
% For example:
 
x>y | x==2 | y==8
 
%% Selecting data using logical statements
%
% We can use logical statements to select data from arrays using the |find|
% command. Type
 
x=[1 3 4 8 6 5 3 2 2]
y=[3 3 1 5 3 2 8 4 2]
condition=y>2
ind=find(condition)
 
%%
% Now the vector |ind| contains the locations of those values of y which
% are greater than 2. Now type 
 
y(1,ind)
 
%%
% MATLAB has returned only those values of y which are greater than 2. Now
% type
 
z=x(1,ind)
 
%%
% We have now stored in the array, z, all the values of x whose 
% corresponding values of y are greater than 2.
%
% Interestingly, the above four lines of code can be combined into a single line
% as follows

z=x(y>2)

%% Using sum, mean, std, max and min
%
% MATLAB has many built in functions of interest. It would be impossible to
% discuss all of them today. But we can look at a few of interest as we go through
% the year.
%
% Enter the following 2D array in the Command Window
 
A=[2 5 8 6; 9 4 3 1; 7 12 10 11]
 
%%
% The |sum| function adds up all of the elements of an array.
%
% |sum(A,1)| sums all the columns of A
 
sum(A,1)
 
%%
% |sum(A,2)| sums all the rows of A
 
sum(A,2)
 
%%
% |sum(A(:),1)| sums everything in A
 
sum(A(:),1)
 
%%
% The |mean| functions works in exactly the same way but finds the
% arithmetic mean, respectively.
 
mean(A,1)
mean(A,2)
mean(A(:),1)
 
%%
% Three other functions also of interest are |max|, |min| and |std|.
%
% |max(A,[],1)| finds the maximum value in each column of A
 
max(A,[],1)
 
%%
% |max(A,[],2)| finds the maximum value in each row of A
 
max(A,[],2)
 
%%
% |max(A(:),[],1)| finds the maximum value in all of A
 
max(A(:),[],1)
 
%%
% Note the inclusion of the |[]|. This is needed as |max| has other
% options, which you are unlikely to consider as you start with MATLAB.
% Therefore we are leaving that argument as empty, i.e., |[]|.
%
% The functions |min| and |std| work in exactly the same way as |max| but
% find the minimum value and the standard deviation, respectively.
%
% To find out more about any function simply type |help| followed by the
% function name in the Command Window. Try typing |help min| in the Command
% Window.
%
%% Generating annual rainfall from monthly rainfall
%
% We will now revisit the Durham weather data set from last week. Recall
% that the data is recorded for every month from 1880 to 2012. The rainfall
% data represents the total depth of rainfall that landed during each
% month. It is not so helpful to look at the statistics of monthly data
% because of the associated effects of seasonality on rain. Instead it
% would be more interesting to study the total rainfall that landed over
% the year. In this exercise we will generate a simple mfile to
% automatically calculate the annual rainfall for each year from 1880 to
% 2012.
%
% Click on "File", "New", "Script". The "Editor" window should now appear.
% Save the file as "MATLABsession2_Assignment.m" in the your "MATLABnotes"
% directory.
%
% At the top of the window type the following:

%This is a script to calculate annual rainfall at Durham from a monthly
%rainfall time-series

%Import Durham weather data from excel
data=readmatrix('DurhamWeather.xlsx','Sheet','Sheet1');

%Extract years and store in yyyy
yyyy=data(:,1);

%Extract months and store in mm
mm=data(:,2);

%Extract monthly rainfall and store in rain
rain=data(:,6);

%Store the first year's worth of data in a separate array and display
rain1880=rain(1:12,1)

%%
% Press the green arrow to run. You should now see the first 12 rainfall
% data in the Command Window. Compare with the excel spreadsheet to verify
% that these are the monthly rainfall data for the year of 1880.
%
% To find the sum for this year, add the following code to your mfile and
% run:

sum(rain1880,1)

%%
% But what if we had to find the sum of every year?
%
% Lets go back to 1880. Add the following code to your mfile and run:

%Choose year of interest
YEAR=1880;
%Obtain the locations of all yyyy which are equal to value of YEAR
condition=yyyy==YEAR;
ind=find(condition);
%Store all rain data for year of interest in rainYEAR
rainYEAR=rain(ind,1);
%Calculate annual rainfall for that year by summation
rainANNUAL=sum(rainYEAR,1)

%%
% We have generalised what we wanted to do for any value stored in |YEAR|.
%
% Now to look at years from 1880 to 2012 we can use a for loop. Add the
% following code to your mfile and run:

%Generate a column vector of the years of interest
YEAR=[1880:2012]';
%Count how many years we are interested in and store in M
M=numel(YEAR);

%Initialise rainANNUAL array before allocating values
rainANNUAL=zeros(M,1);
%Perform a loop for each year of interest
for i=1:M
    %Obtain the locations of all yyyy which are equal to the ith YEAR
    condition=yyyy==YEAR(i,1);
    ind=find(condition);
    %Store all rain data for year of interest in rain YEAR
    rainYEAR=rain(ind,1);
    %Calculate annual rainfall for that year by summation and store in the
    %ith value of rainANNUAL
    rainANNUAL(i,1)=sum(rainYEAR,1);
end

%%
% Interestingly, the code would work just as well without the
% |rainANNUAL=zeros(M,1)| step. However, MATLAB
% will compute the process faster if you specify the size of an array
% before allocating values in loops. The |zeros(M,1)| command simply gives
% an array of zeros with M rows and 1 column.

%% Classroom assignment
%
% 1) Extend your code further by finding the a) maximum, b) minimum, c) mean 
% and d) standard deviation of the annual rainfall in Durham for the 
% period 1980 to 2011, respectively.
%
% 2) Extend the code to count how many years during that period had
% above average rainfall.
%
% 3) Use MATLAB to generate an appropriately labelled bar chart
% displaying the annual rainfall for the period 1980 to 2011. Type |help
% bar| to find out how to do this. Hint: it is very similar to using the
% |plot| command you used last week. Use the bar chart to check your ansers
% to the previous two tasks.
%
% An example MATLAB code, containing all of the instructions requested
% above, is given in <https://github.com/samathiasuk/MATLABnotes/blob/main/MATLABsession2_Assignment.m MATLABsession2_Assignment.m>.
