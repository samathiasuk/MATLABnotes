function MATLABsession3_Assignment
%This file contains instructions requested by MATLAB Notes - Session 3

%Import Durham weather data from excel
data=readmatrix('DurhamWeather.xlsx','Sheet','Sheet1');
%Extract years and store in yyyy
yyyy=data(:,1);
%Extract months and store in mm
mm=data(:,2);
%Extract monthly rainfall and store in rain
rain=data(:,6);
%Extract monthly mean daily maximum temperature (deg C) and store in tmax
tmax=data(:,3);
%Generate a column vector of the years of interest
YEAR=[1880:2011]';
%Apply sub function to get CDF for annual mean monthly rainfall
[RankedRain,P]=GetCDF(YEAR,yyyy,rain);
%Apply sub function to get CDF for annual mean daily max temperature
[RankedTmax,P]=GetCDF(YEAR,yyyy,tmax);
%Get thoretical normal CDF for rain
[xRain,PnormalRain]=GetNormalDist(RankedRain);
%Get thoretical normal CDF for temp
[xTmax,PnormalTmax]=GetNormalDist(RankedTmax);
figure(1)
clf
subplot(1,2,1)
%Plot data as x-y scatter
plot(RankedRain,P,'.',xRain,PnormalRain)
xlabel('Annual mean rainfall (mm/month)')
ylabel('Probability of non-exceedance')
title('Rainfall in Durham from 1880 to 2012')
legend('Empirical','Normal')
subplot(1,2,2)
plot(RankedTmax,P,'.',xTmax,PnormalTmax)
xlabel('Annual mean daily max temp (^oC)')
ylabel('Probability of non-exceedance')
title('Temperature in Durham from 1880 to 2012')
legend('Empirical','Normal')

function [RankedRain,P]=GetCDF(YEAR,yyyy,rain)
%Count how many years we are interested in and store in M
M=numel(YEAR);
%Perform a loop for each year of interest
for i=1:M
    %Obtain the locations of all yyyy which are equal to the ith YEAR
    condition=yyyy==YEAR(i,1);
    ind=find(condition);
    %Store all rain data for year of interest in rainYEAR
    rainYEAR=rain(ind,1);
    %Calculate annual mean rainfall for that year
    %and store as the ith value of rainANNUAL
    rainANNUAL(i,1)=mean(rainYEAR,1);
end
%Rank the rain from smallest to largest
RankedRain=sort(rainANNUAL,1);
%Count the number of data points
M=numel(RankedRain);
%Generate the rank numbers (i.e., numbers from 1 to M in increments of 1)
m=[1:M]';
%Calculate probability of non-exceedance using Weibull plotting position
P=m/(M+1);
function [x,Pnormal]=GetNormalDist(RankedRain)
%Calculate mean
mu=mean(RankedRain,1);
%Calculate a standard deviation
sigma=std(RankedRain,[],1);
%Generate a set of 50 equally spaced points from the min to the max
rainMIN=min(RankedRain,[],1);
rainMAX=max(RankedRain,[],1);
%Type "help linspace" in Command Window to find out more.
x=linspace(rainMIN,rainMAX,50)';
%Evaluate Normal CDF
u=(mu-x)/sigma/sqrt(2);
Pnormal=erfc(u)/2;