%% Initialization

clear all
close all
clc

% This script will analyze the data from the EQ Qualtrics survey for
% Daniel's thesis.
% Daniel Sazhin
% 03/26/2018

[~,~,data] = xlsread('Thesis_Data_03262.csv') 

%% Sorting

rawdata = data(4:end,:) % Take all numerical calues
rawdata(:,1:17) = []
