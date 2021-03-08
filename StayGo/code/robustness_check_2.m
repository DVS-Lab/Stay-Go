%% Initialization

clear all;
close all;
clc;


%% Briefing

% This script analyzes Stay/Go data. It is organized by the hypotheses
% assessed. See each subsection for more detailed information.

% Daniel Sazhin
% 07-30-2020
% DVS Lab
% Temple University

load Survey_Data_robust.mat 

N=121;

%% Hypotheses

% Given exponentially varying information (defined as a steeper/shallower curve), 
% how well do participants predict a future trend (i.e.,the next direction that the curve will take)?
%
% In this study, we expect that if a trial is trending better, with a participant 
% expecting to make more money, that their estimate of the final value of the
% trend should be a factor in their choice. Secondly, the extent to which the 
% information is exponential, represented by a growth factor will affect a 
% participant’s choice. Depending on the final value and growth factor of a 
% given trial, the participant will decide whether they will earn more by 
% staying through all ten turns of a trial or leave at an earlier turn.
%
% H1: Participants will behave suboptimally; they will make less money than 
% chance. (Optimality is defined as randomly selecting a strategy of staying 
% all the way through or quitting on any given trial.)
%
% H2: Higher growth factors will lead to subjects waiting longer.
%
% H2.1: Higher growth factors lead to participants waiting more turns before quitting.
% H2.2: Higher growth factors lead to participants quitting when they should stay.
% H2.3: Higher final values lead to participants more likely to stay than to quit.
% H2.4: There will be an interaction effect between final values and growth 
% factors, with higher growth factors minimizing the likelihood of subjects staying through the end.
%
% H3: Higher growth factors will lead to subjects making less money.
%
% H4: Individual difference measures will be associated with behavior
%
% H4.1: Higher risk aversion will be associated with participants making less money.
% H4.2: Higher loss aversion will be associated with participants making less money.
% H4.3: Higher intolerance of uncertainty will be associated with participants making less money.
% H4.4: Greater mood symptoms will be associated with participants making less money.
% H4.5: Greater substance use will be associated with participants making less money.


%% H1: Are subjects optimal?

% We defined suboptimality as making less than always staying until 10.

% Find average earnings for each participant

earnings= [];

for ii = 1:N
   
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
participant_earnings =(Participant(:,4));
earnings = [participant_earnings; earnings];
average_earnings = mean(earnings);
    
end

[H,P,CI,stats] = ttest(average_earnings,10); % performs a t-test of the hypothesis that the data in X come from a distribution with mean M.  M must be a scalar.

figure
h = histogram(average_earnings);
counts = h.Values;
h.NumBins = 11;
ax = gca;
ax.FontSize = 9;
xlabel ('Average Trial Earnings ($)','FontSize', 16)
ylabel ('Frequency','FontSize', 16)
set(gca,'box','off')
set(gcf,'color','w');


%% Other strategies 

[n,t,rawdata] = xlsread('Exponentials.csv');
testdata = cell2mat(rawdata);

% Always go to 10

EV_AlwaysTen = mean(testdata(:,10));

% Always leave at 1 

EV_AlwaysOne = 10 - mean(testdata(:,1));

% Randomly leave at move 1-10

% Let's do a monte carlo simulation

total_simulated_earnings = [];

for yyy = 1:1000 % number of simulations
    subject_earnings  = [];
    for xxx = 1:length(testdata) % go through the triials
        r = randi([1 10]); % randomly select a turn to leave
        if r < 10
            turn_earnings = 10 - testdata(xxx,r); % earnings
            subject_earnings = [subject_earnings; turn_earnings];
        end
        if r == 10
            turn_earnings = testdata(xxx,r); % earnings
            subject_earnings = [subject_earnings; turn_earnings];
        end
        
    end
    
    average_earnings_sim = mean(subject_earnings);
    total_simulated_earnings = [total_simulated_earnings; average_earnings_sim];
end
    
mean_total_simulated_earnings = mean(total_simulated_earnings);
figure
h = histogram(total_simulated_earnings(:));
counts = h.Values;
h.NumBins = 11;
ax = gca;
ax.FontSize = 9;
xlabel ('total simulated earnings ($) (random)','FontSize', 16)
ylabel ('Frequency','FontSize', 16)
set(gca,'box','off')
set(gcf,'color','w');

% Assumption: A computer could fit an exponetial function almost perfectly at the 4th turn
% 50 percent of the time, leaving at the fourth turn is $9.6745.
% 50 percent of the time, predicting the final value as higher is $12.5
% EV of optimal strategy is $11.09

%% H4: No association between earnings and scales.

% [R,P] = corrcoef(average_earnings,Uncertainty_data_final);
% [R,P] = corrcoef(average_earnings,AADIS_data_final);
% [R,P] = corrcoef(average_earnings,SevenUp);
% [R,P] = corrcoef(average_earnings,SevenDown);


%% H2: Higher growth factors will lead to subjects waiting longer.

Growth_Total = [];
Final_Value_Total = [];
Turn_Left_Total = [];
Earnings_Total = [];

betas_first = [];
stats_save = [];

for ii = 1:N
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

x1 = reshape(Final_Value,[],1);
x1= zscore(x1);
x2 = reshape(Growth,[],1);
x2= zscore(x2);
y = reshape(Turn_Left,[],1);
y = zscore(y);

X = [ones(size(x1)) x1 x2 x1.*x2];
[b,bint,r,rint,stats] = regress(y,X);    % Removes NaN data

betas_first = [betas_first,b];
stats_save = [stats_save;stats];

end

[A,B] = size(betas_first);
save_ttest = [];
for jj = 1:A
    row = betas_first(jj,:);
    [H,P,CI,stats] = ttest(row,0);
    save_ttest = [save_ttest, P];
end

B2Er = std(betas_first(2,:)') / sqrt(length(betas_first(2,:)'));
B3Er = std(betas_first(3,:)') / sqrt(length(betas_first(3,:)'));
B4Er = std(betas_first(4,:)') / sqrt(length(betas_first(4,:)'));

err = [B2Er,B3Er,B4Er] * 2;

data = [mean(betas_first(2,:));  mean(betas_first(3,:)); mean(betas_first(4,:))]';
x = linspace(1,3,3);
figure
bar(x,data)    
ax = gca;
ax.FontSize = 12;
xlim ([.5 3.5]);
box off
xlabel ('Regressors', 'FontSize', 16);
ylabel  ('Z-Standardized Beta Weight', 'FontSize', 16);
set(gcf,'color','w');
set(gca, 'XTick', 1:3, 'XTickLabels', {'Final Value','Growth Factor','Interaction'})
title('Effect of Growth Factors and Final Values on Turn Left')

hold on

% Standard Error 

er = errorbar(x,data,err); %errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
hold off

saveas(gcf,'Bar_Left.png')

% Find subject betas 

[Betas_min,Index__min] = min(betas_first');
[Betas_max,Index__max] = max(betas_first');

% 2 = final value
% 3 = growth factor
% 4 = interaction

%% let's look at a subject with low and high final values.

low_beta_final = Index__min(2);

ii = low_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)
ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('low beta participant for final value (H2 Model)')
ylim([5 15])

high_beta_final = Index__max(2);

ii = high_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)
ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('high beta participant for final value (H2 Model)')
ylim([5 15])


%% let's look at a subject with low and high growth values.

low_beta_final = Index__min(3);

ii = low_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)
ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('low beta participant for growth factor (H2 Model)')
ylim([5 15])


high_beta_final = Index__max(3);

ii = high_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)
ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('high beta participant for growth factor (H2 Model)')
ylim([5 15])


%% H3: Higher growth factors will lead to subjects making less money.

Growth_Total = [];
Final_Value_Total = [];
Turn_Left_Total = [];
Earnings_Total = [];

betas = [];
stats_save = [];

for ii = 1:N;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

x1 = reshape(Final_Value,[],1);
x1 = zscore(x1);
x2 = reshape(Growth,[],1);
x2 = zscore(x2);
y = reshape(Earnings,[],1);
y = zscore(y);

X = [ones(size(x1)) x1 x2 x1.*x2];
[b,bint,r,rint,stats] = regress(y,X);    % Removes NaN data

betas = [betas,b];
stats_save = [stats_save;stats];

end

[A,B] = size(betas);
save_ttest = [];
for jj = 1:A
    row = betas(jj,:);
    [H,P,CI,stats] = ttest(row,0);
    save_ttest = [save_ttest, P];
end

B2Er = std(betas(2,:)') / sqrt(length(betas(2,:)'));
B3Er = std(betas(3,:)') / sqrt(length(betas(3,:)'));
B4Er = std(betas(4,:)') / sqrt(length(betas(4,:)'));

err = [B2Er,B3Er,B4Er] * 2;

data = [mean(betas(2,:));  mean(betas(3,:)); mean(betas(4,:))]';
x = linspace(1,3,3);
figure
bar(x,data)    
ax = gca;
ax.FontSize = 12;
xlim ([.5 3.5]);
box off
xlabel ('Regressors', 'FontSize', 16);
ylabel  ('Z-Standardized Beta Weights', 'FontSize', 16);
set(gcf,'color','w');
set(gca, 'XTick', 1:3, 'XTickLabels', {'Final Value','Growth Factor','Interaction'})
title('Effect of Growth Factors and Final Values on Earnings')

hold on

% Standard Error 

er = errorbar(x,data,err); %errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1;
hold off

saveas(gcf,'Bar_Earnings.png')

%% let's look at a subject with low and high final values.

low_beta_final = Index__min(2);

ii = low_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)
ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('low beta participant for final value ((H3 Model)')
ylim([5 15])

high_beta_final = Index__max(2);

ii = high_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)

ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('high beta participant for final value (H3 Model)')
ylim([5 15])

%% let's look at a subject with low and high interaction factors.

low_beta_final = Index__min(4);

ii = low_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)

ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('low beta participant for interaction factor (H3 Model)')
ylim([5 15])

high_beta_final = Index__max(4);

ii = high_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Turn_Left,Earnings)

ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('high beta participant for interaction factor (H3 Model)')
ylim([5 15])


