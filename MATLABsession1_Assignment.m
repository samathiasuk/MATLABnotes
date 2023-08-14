%This file contains instructions requested by MATLAB Notes - Session 1
 
%Just to get started, make a row vector of numbers
A=[4 6 8 5 1 3 2]
%Multiply A by 2 and store as B
B=2*A

%Import Durham weather data from Excel
data=readmatrix('DurhamWeather.xlsx','Sheet','Sheet1');

%Determine the size of the data
size(data)

%Extract years and store in yyyy
yyyy=data(:,1);

%Extract months and store in mm
mm=data(:,2);

%Extract monthly mean daily minimum temperature and store in tmin
tmin=data(:,4);

%Extract monthly number of days with air-frost and store in af
af=data(:,5);

%Combine yyyy and mm to form a decimalised year
tYr=yyyy+mm/12;

%Make a figure
figure(1)
%Clear figure(1) if it alreay exists
clf
%Plot time against minimum temperature and days with air-frost
plot(tYr,tmin,tYr,af)

%Provide a legend, axes labels and a title
legend('Minimum temperature in ^oC','Monthly days of air-frost')
xlabel('Year')
ylabel('See legend')
title('The weather in Durham')

%Limit the axes to a range of interest
axis([2000 2013 -5 25])

%Classroom assignment

% This script imports and plots permeability data and porosity data from a range of
% North Sea sandstones and shales.

%Import Sandstone data
SandstoneData=readmatrix('NorthSeaRocks.xlsx','Sheet','Sandstones');

%Import Shale data
ShaleData=readmatrix('NorthSeaRocks.xlsx','Sheet','Shales');

%Extract sandstone porosities and permeabilities 
SandstonePoro=SandstoneData(:,1:2:end);
SandstonePerm=SandstoneData(:,2:2:end);

%Extract shale porosities and permeabilities 
ShalePoro=ShaleData(:,1:2:end);
ShalePerm=ShaleData(:,2:2:end);

%Use figure(2) to make a new figure
figure(2) 
%Clear figure(1) if it alreay exists
clf
%Plot data on double log axes
loglog(SandstonePoro,SandstonePerm,'o',ShalePoro,ShalePerm,'x')
xlabel('Porosity [frac]')
ylabel('Permeability [mD]')
legend('Bunter Sandstone','Rotliegend Sandstone','Carboniferous Sandstone',...
    'Jurassic Shale','Rotliegend Shale','Carboniferous Shale','location','northwest')
