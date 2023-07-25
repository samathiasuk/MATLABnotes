%This file contains instructions requested by MATLAB Notes - Session 2

%The objective is to calculate annual rainfall at Durham from a monthly 
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
%Find the sum
sum(rain1880,1)

%Choose year of interest
YEAR=1880;
%Obtain the locations of all yyyy which are equal to value of YEAR
condition=yyyy==YEAR;
ind=find(condition);
%Store all rain data for year of interest in rainYEAR
rainYEAR=rain(ind,1);
%Calculate annual rainfall for that year by summation
rainANNUAL=sum(rainYEAR,1)
    
%Generate a column vector of the years of interest
YEAR=[1980:2012]';
%Count how many years we are interested in and store in M
M=numel(YEAR);
%Perform a loop for each year of interest
for i=1:M
    %Obtain the locations of all yyyy which are equal to the ith YEAR
    condition=yyyy==YEAR(i,1);
    ind=find(condition);
    %Store all rain data for year of interest in rainYEAR
    rainYEAR=rain(ind,1);
    %Calculate annual rainfall for that year by summation
    %and store as the ith value of rainANNUAL
    rainANNUAL(i,1)=sum(rainYEAR,1);
end

%Classroom asignment:

%For the period of 1980 to 2011 in Durham calculate the:
%1a) Maximum annual rainfall
MaxRain=max(rainANNUAL,[],1)
%1b) Minimum annual rainfall
MinRain=min(rainANNUAL,[],1)
%1c) Mean annual rainfall
MeanRain=mean(rainANNUAL,1)
%1d) Standard deviation of the annual rainfall
STDofRain=std(rainANNUAL,[],1)

%2) Calculate number of years where rainfall is above average
%Find the locations of the above average years
condition=rainANNUAL>MeanRain;
ind=find(condition);
%Use numel to count how many there are
NumberOfAboveAverageYears=numel(ind)

%3) Plot the results as a barchart
figure(1)
bar(YEAR,rainANNUAL)
axis([1980 2012 0 1000])
xlabel('Year')
ylabel('Annual rainfall (mm)')
title('Annual rainfall at Durham from 1980 to 2011')
