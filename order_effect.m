clear all
close all
clc

[n,t,rawdata] = xlsread('Thesis_Data_0403_leave.csv'); % This converts the CSV file into a matlab file. 

% Find FL_163_DO

temp = find(strcmp('FL_163_DO',rawdata(1,:)));
FL_163_DO = rawdata(:,temp);
order =(FL_163_DO(4:end,:)); % This is the subset of data we will use in matrix form
save('Order','order','FL_163_DO')