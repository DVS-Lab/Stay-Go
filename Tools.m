clear all
close all
clc

%% Useful tools for the future!

%% Bar Plot

load Behavioral_Data.mat

% Calculate Error Bars by taking two Standard Errors.

% Standard Error for error bars

% ErrBar1 = std(..)/ sqrt(length(..))
%err = [ErrBar1,etc...] * 2; % This squares the error and saves the error
%vector.


figure
% x = linspace(0,.6,7) % Set linspace to the X axis. This sets the X axis
% of the bar chart.
% y = [...]; Y is the data.
%errhigh = [ErrBar1 + mean(first_bar)... etc];
%errlow  = [ErrBar1 - mean(first_bar)... etc];
figure

% 
bar(x,y)    
ax = gca
ax.FontSize = 12 % Sets fontsize of the axis
xlim ([-.05 .65]) % Sets x limits
box off % Turns off the ticks
xlabel ('Offers', 'FontSize', 16); % X label with font
ylabel  ('P(Accept)', 'FontSize', 16); % Y label with font
set(gcf,'color','w'); % Makes background white instead of gray.

hold on

er = errorbar(x,y,err) % Error bar inserted
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

saveas(gcf,'Fig4.png') % Saves the figure.

%% Scatterplot

% [R,P] = corrcoef(X, Y)
figure
scatter(X, Y, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5) % These settings change the dots.
ax = gca
ax.FontSize = 12 % Font size of the units.
xlabel ('X', 'FontSize', 16); % Labels and change their fontsizes.
ylabel  ('Y', 'FontSize', 16);
lsline % Makes a best fit line.
set(gcf,'color','w');

saveas(gcf,'Fig10.png')

%% Histogram 

figure
h = histogram(X);
counts = h.Values;
h.NumBins = 12
ax = gca
ax.FontSize = 12
xlabel ('X','FontSize', 16)
ylabel ('Frequency','FontSize', 16)
set(gca,'box','off')
set(gcf,'color','w');

saveas(gcf,'Fig1.png')

%% Logistic Regression

[B,DEV,STATS] = mnrfit(Data1, Data2)

% Open stats for interpretation. Use t and p values.

%% Correlation

[R,P] = corrcoef(Data1, Dat2)

%% T-test

[H,P_UG_A_Earnings,CI,STATS] = ttest2(Data1(:,1), Data2(:,1));

