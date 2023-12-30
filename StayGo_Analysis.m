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

load Survey_Data.mat

[N,M] = size(Uncertainty_data_final); % N is the number of subjects

%% H1: Are subjects optimal?

% We defined suboptimality as making less than always staying until 10.

% Find average earnings for each participant

earnings= [];

for ii = 1:N
    filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
    
    Participant = csvread(filename,1,0);
    participant_earnings =(Participant(:,4));
    earnings = [participant_earnings, earnings];
    average_earnings = mean(earnings)';
    
end

[H,P,CI,stats] = ttest(average_earnings,10); % performs a t-test of the hypothesis that the data in X come from a distribution with mean M.  M must be a scalar.

figure
h = histogram(average_earnings(:));
counts = h.Values;
h.NumBins = 11;
ax = gca;
ax.FontSize = 9;
xlabel ('Average Trial Earnings ($)','FontSize', 16)
ylabel ('Frequency','FontSize', 16)
set(gca,'box','off')
set(gcf,'color','w');
saveas(gcf,'H1.svg')


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


%% Manic-Depressive Symptoms

figure
scatter(SevenUp, average_earnings, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca;
ax.FontSize = 12
xlabel ('Mood Symptoms (7Up)', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
i = lsline;
i.LineWidth = 3;
i.LineStyle = '--'
i.Color = [1 0 0];
set(gcf,'color','w');

%% H4: No association between earnings and scales.

[R,P] = corrcoef(average_earnings,Uncertainty_data_final);

figure
scatter(average_earnings,Uncertainty_data_final, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Earnings', 'FontSize', 16);
ylabel  ('Intolerance of Uncertainty', 'FontSize', 16);
i = lsline;
i.LineWidth = 3;
i.LineStyle = '--'
i.Color = [1 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_H1_IUD.svg')

[R,P] = corrcoef(average_earnings,AADIS_data_final);

figure
scatter(average_earnings,AADIS_data_final, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca;
ax.FontSize = 12;
xlabel ('Earnings', 'FontSize', 16);
ylabel  ('Substance Use (AADIS)', 'FontSize', 16);
i = lsline;
i.LineWidth = 3;
i.LineStyle = '--';
i.Color = [1 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_H1_AADIS.svg')

[R,P] = corrcoef(average_earnings,SevenUp);

figure
scatter(average_earnings,SevenUp, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Earnings', 'FontSize', 16);
ylabel  ('Mood Symptoms (7Up)', 'FontSize', 16);
i = lsline;
i.LineWidth = 3;
i.LineStyle = '--'
i.Color = [1 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_H1_SevenUp.svg')

[R,P] = corrcoef(SevenDown, average_earnings);

figure
scatter(average_earnings, SevenUp, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Earnings', 'FontSize', 16);
ylabel  ('Mood Symptoms (7Up)', 'FontSize', 16);
i = lsline;
i.LineWidth = 3;
i.LineStyle = '-'
i.Color = [1 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_H1_SevenDown.svg')

[R,P] = corrcoef(average_earnings,save_risk);

figure
scatter(average_earnings,save_risk, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Earnings', 'FontSize', 16);
ylabel  ('Risk', 'FontSize', 16);
i = lsline;
i.LineWidth = 3;
i.LineStyle = '--'
i.Color = [1 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_H1_Risk.svg')

[R,P] = corrcoef(average_earnings,save_loss);

figure
scatter(average_earnings,save_loss, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Earnings', 'FontSize', 16);
ylabel  ('Risk', 'FontSize', 16);
i = lsline;
i.LineWidth = 3;
i.LineStyle = '--'
i.Color = [1 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_H1_Loss.svg')


% [R,P] = corrcoef(average_earnings,Risk_Indifference);
% [R,P] = corrcoef(average_earnings,Loss_Indifference);


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

err = [B2Er,B3Er,B4Er] * 2;

er = errorbar(x,data,err); %errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
hold off

saveas(gcf,'Bar_Left.svg')

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
scatter(Turn_Left,Growth)
ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Growth', 'FontSize', 16);
set(gcf,'color','w');
title('low beta participant for final value (H2 Model)')

x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);

figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('low beta participant for final value (H2 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);


high_beta_final = Index__max(2);

ii = high_beta_final;
filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
Participant = csvread(filename,1,0);
Growth = Participant(:,1);
Final_Value = Participant(:,2);
Turn_Left = Participant(:,3);
Earnings = Participant(:,4);

figure
scatter(Growth,Earnings)
ax = gca;
ax.FontSize = 12;
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
set(gcf,'color','w');
title('high beta participant for final value (H2 Model)')
ylim([5 15])

x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);

figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('high beta participant for final value (H2 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);


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

x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);

figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('low beta participant for growth (H2 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);




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


x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);


figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('high beta participant for growth (H2 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);

%% Exploratory: analysis of IDs with Effect of regressors on turn left model?

% Is there an association of final values regressor in the (effect of regressors on earnings) model.

Final_values_beta =  betas_first(2,:)'; % Betas associated with Final Values

[R,P] = corrcoef(Final_values_beta,SevenDown); % No
[R,P] = corrcoef(Final_values_beta,SevenUp); % Yes!

figure
scatter(Final_values_beta,SevenDown, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Final Values', 'FontSize', 16);
ylabel  ('Seven Down', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_Exploratory_Sevendown.svg')

figure
scatter(Final_values_beta,SevenUp, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Final Values', 'FontSize', 16);
ylabel  ('Seven Up', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_Exploratory_Sevenup.svg')

[R,P] = corrcoef(Final_values_beta,AADIS_data_final); % Nein
[R,P] = corrcoef(Final_values_beta,Uncertainty_data_final); % Yes!

figure
scatter(Final_values_beta,Uncertainty_data_final, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Final Values', 'FontSize', 16);
ylabel  ('Uncertainty', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

%[R,P] = corrcoef(Final_values_beta,Risk_Indifference); % No
% [R,P] = corrcoef(Final_values_beta,Loss_Indifference);

% Is there an association with Growth Factor in Turn Left with IDs?

Growth_values_beta =  betas_first(3,:)'; % Betas associated with Final Values

[R,P] = corrcoef(Growth_values_beta,SevenDown); % No
[R,P] = corrcoef(Growth_values_beta,SevenUp); % No
[R,P] = corrcoef(Growth_values_beta,AADIS_data_final); % Nein
[R,P] = corrcoef(Growth_values_beta,Uncertainty_data_final); % Nyet
%[R,P] = corrcoef(Growth_values_beta,Risk_Indifference); % No
%[R,P] = corrcoef(Growth_values_beta,Loss_Indifference); % No, but with
% more data is promising.

% Is there an association with Interaction Factor in Earnings with IDs?

Interaction_values_beta =  betas_first(4,:)'; % Betas associated with Final Values

[R,P] = corrcoef(Interaction_values_beta,SevenDown); % No
[R,P] = corrcoef(Interaction_values_beta,SevenUp); % No
[R,P] = corrcoef(Interaction_values_beta,AADIS_data_final); % Nein
[R,P] = corrcoef(Interaction_values_beta,Uncertainty_data_final); % Nyet
% [R,P] = corrcoef(Interaction_values_beta,Risk_Indifference); % No
% [R,P] = corrcoef(Interaction_values_beta,Loss_Indifference); % No


%% Apply multi-level model for H2

% We will run a multilevel model, wherein the dependent variable is the number 
% of turns a participant chooses to stay (turns) as a function of the final value (value), 
% growth factor (growth) and interaction between these two factors.
% We will include random effects of participant, and random slopes for final
% value and growth nested within participant 
% (in R notation): turns ~ value * growth + (value | ID) + (growth | ID).
% 
Growth = [];
Final_Value = [];
Turn_Left = [];
Earnings = [];
ID = [];

for ii = 1:N
    filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
    Participant = csvread(filename,1,0);
    
    Growth_Participant = Participant(:,1);
    Final_Value_Participant = Participant(:,2);
    Turn_Left_Participant = Participant(:,3);
    Earnings_Participant = Participant(:,4);
    
    [A,B] = size(Participant);
    Subject = ii*ones(A,1);
    
    Growth = [Growth; Growth_Participant];
    Final_Value = [Final_Value; Final_Value_Participant];
    Turn_Left = [Turn_Left; Turn_Left_Participant];
    Earnings = [Earnings; Earnings_Participant];
    ID = [ID; Subject];
    
end
    
tb1 = array2table(zscore(Growth),'VariableNames', {'growth'}); 
tb2 = array2table(zscore(Turn_Left),'VariableNames', {'turns'});
tb3 = array2table(zscore(Final_Value),'VariableNames', {'value'});
tb4 = array2table(ID,'VariableNames', {'ID'});

%tb1 = zscore(tb1);
%tb2 = zscore(tb2);
%tb3 = zscore(tb3);
%tb4 = zscore(tb4);
mixed_effects = [tb2, tb3, tb1, tb4];  
fixed_effects = [tb2, tb3, tb1]; 

random_slope_model = fitlme(mixed_effects,'turns ~ value * growth + (value | ID) + (growth | ID)')
null_effect_model = fitlme(mixed_effects,'turns ~ 1 + (value | ID) + (growth | ID)')

intercept_results = random_slope_model.Coefficients(1,:);
value_results = random_slope_model.Coefficients(2,:);
growth_results = random_slope_model.Coefficients(3,:);
interaction_results = random_slope_model.Coefficients(4,:);

data = [value_results.Estimate;  growth_results.Estimate; interaction_results.Estimate]';
err = [value_results.SE,growth_results.SE,interaction_results.SE] * 2;
x = linspace(1,3,3);
figure
bar(x,data)
ax = gca;
ax.FontSize = 12;
xlim ([.5 3.5]);
box off
xlabel ('Regressors', 'FontSize', 16);
ylabel  ('Beta Estimate', 'FontSize', 16);
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

saveas(gcf,'Bar_Left.svg')


compare(null_effect_model,random_slope_model);
[~,~,RE] = randomEffects(random_slope_model);
[~,~,stats] = covarianceParameters(random_slope_model);

[r, LB, UB, F, df1, df2, p] = ICC([zscore(Growth),zscore(Final_Value),zscore(Turn_Left)],'C-1', .05, 0)

testmat= [zscore(Growth),zscore(Final_Value),zscore(Turn_Left)]
collintest(mixed_effects,Plot="on");

figure
collintest(testmat,Plot="on", ...
    TolIdx=10,TolProp=0.5);

[T,P,df] = BPtest(testmat)


% apply VIF test to all the independent variables
test_vars = [zscore(Growth),zscore(Final_Value)];
interaction = (test_vars(:,1) - mean(test_vars(:,1))) .* (test_vars(:,2)) - mean(test_vars(:,2));
independent_vars = [test_vars, interaction];

test = vif(independent_vars)
%% H3: Higher growth factors will lead to subjects making less money.

Growth_Total = [];
Final_Value_Total = [];
Turn_Left_Total = [];
Earnings_Total = [];

betas = [];
betas_first = [];
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
    % x3 = reshape(Turn_Left,[],1);
    % x3 = zscore(x3);
    y = reshape(Earnings,[],1);
    y = zscore(y);
    
    X = [ones(size(x1)) x1 x2 x1.*x2];
    % X = [ones(size(x1)) x1 x2 x1.*x2 x2];
    [b,bint,r,rint,stats] = regress(y,X);    % Removes NaN data
    
    betas = [betas,b];
    betas_first = [betas_first,b];
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

saveas(gcf,'Bar_Earnings.svg')

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

x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);

figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('low beta participant for final value (H3 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);

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

x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);

figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('high beta participant for final value (H2 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);
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

x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);

figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('low beta participant for interaction (H3 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);

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

x = Turn_Left';
y = Earnings';
z = Growth';
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 10)';
yv = linspace(min(y), max(y), 10)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);

figure
scatter3(x,y,z, 'filled')
hold on
mesh(X, Y, Z, 'FaceAlpha', 0.5)
hold off
view(-120, 35)
title(sprintf('high beta participant for interaction (H3 Model)', B))
xlabel ('Turn Left', 'FontSize', 16);
ylabel  ('Earnings', 'FontSize', 16);
zlabel ('Growth', 'Fontsize', 16);


%% Apply multi-level model for H3

% We will run a multilevel model, wherein the dependent variable is participant earnings
% as a function of the final value (value), 
% growth factor (growth) and interaction between these two factors.
% We will include random effects of participant, and random slopes for final
% value and growth nested within participant 
% (in R notation): turns ~ value * growth + (value | ID) + (growth | ID).
% 

Growth = [];
Final_Value = [];
Turn_Left = [];
Earnings = [];
ID = [];

for ii = 1:N
    filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
    Participant = csvread(filename,1,0);
    
    Growth_Participant = Participant(:,1);
    Final_Value_Participant = Participant(:,2);
    Turn_Left_Participant = Participant(:,3);
    Earnings_Participant = Participant(:,4);
    
    [A,B] = size(Participant);
    Subject = ii*ones(A,1);
    
    Growth = [Growth; Growth_Participant];
    Final_Value = [Final_Value; Final_Value_Participant];
    Turn_Left = [Turn_Left; Turn_Left_Participant];
    Earnings = [Earnings; Earnings_Participant];
    ID = [ID; Subject];
    
end


tb1 = array2table(zscore(Growth),'VariableNames', {'growth'}); 
tb2 = array2table(zscore(Earnings),'VariableNames', {'earnings'});
tb3 = array2table(zscore(Final_Value),'VariableNames', {'value'});
tb4 = array2table(ID,'VariableNames', {'ID'});

% [N,~] = size(tb1);
% ones_sample = ones(N,1);

mixed_effects = [tb2, tb3, tb1, tb4];  
fixed_effects = [tb2, tb3, tb1];  
random_slope_model = fitlme(mixed_effects,'earnings ~ value * growth + (value | ID) + (growth | ID)')
null_model = fitlme(mixed_effects,'earnings ~ 1 + (value | ID) + (growth | ID)');

intercept_results = random_slope_model.Coefficients(1,:);
value_results = random_slope_model.Coefficients(2,:);
growth_results = random_slope_model.Coefficients(3,:);
interaction_results = random_slope_model.Coefficients(4,:);

data = [value_results.Estimate;  growth_results.Estimate; interaction_results.Estimate]';
err = [value_results.SE,growth_results.SE,interaction_results.SE] * 2;
x = linspace(1,3,3);
figure
bar(x,data)
ax = gca;
ax.FontSize = 12;
xlim ([.5 3.5]);
box off
xlabel ('Regressors', 'FontSize', 16);
ylabel  ('Beta Estimate', 'FontSize', 16);
set(gcf,'color','w');
set(gca, 'XTick', 1:3, 'XTickLabels', {'Final Value','Growth Factor','Interaction'})
title('Effect of Growth Factors and Final Values on Earnings (Multilevel)')

hold on

% Standard Error

er = errorbar(x,data,err); %errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
hold off

saveas(gcf,'Bar_Earnings.svg')


% Z-scoring changes the effect.


compare(null_model,random_slope_model);
[~,~,RE] = randomEffects(random_slope_model);
[~,~,stats] = covarianceParameters(random_slope_model);

[r, LB, UB, F, df1, df2, p] = ICC([zscore(Growth),zscore(Final_Value),zscore(Earnings)],'C-1', .05, 0)

testmat= [zscore(Growth),zscore(Final_Value),zscore(Earnings)];
collintest(mixed_effects,Plot="on");

figure
collintest(testmat,Plot="on", ...
    TolIdx=1,TolProp=0.5);


% apply VIF test to all the independent variables
test_vars = [zscore(Growth),zscore(Final_Value)];
interaction = (test_vars(:,1) - mean(test_vars(:,1))) .* (test_vars(:,2)) - mean(test_vars(:,2));
independent_vars = [test_vars, interaction];

test = vif(independent_vars)

%% Exploratory: analysis of IDs with Effect of regressors on earnings model?

% Is there an association of final values regressor in the (effect of regressors on earnings) model.

Final_values_beta =  betas_first(2,:)'; % Betas associated with Final Values

[R,P] = corrcoef(Final_values_beta,SevenDown); % No
[R,P] = corrcoef(Final_values_beta,SevenUp); % No...
[R,P] = corrcoef(Final_values_beta,AADIS_data_final); % Nein
[R,P] = corrcoef(Final_values_beta,Uncertainty_data_final); % No!
[R,P] = corrcoef(Final_values_beta,AADIS_data_final); % No.
[R,P] = corrcoef(Final_values_beta,save_risk); % No
[R,P] = corrcoef(Final_values_beta,save_loss); % No

% Is there an association with Growth Factor in Turn Left with IDs?

Growth_values_beta =  betas_first(3,:)'; % Betas associated with Final Values

[R,P] = corrcoef(Growth_values_beta,SevenDown); % No
[R,P] = corrcoef(Growth_values_beta,SevenUp); % No
[R,P] = corrcoef(Growth_values_beta,AADIS_data_final); % Nein
[R,P] = corrcoef(Growth_values_beta,Uncertainty_data_final); % Yes! 
%[R,P] = corrcoef(Growth_values_beta,Risk_Indifference); % No
%[R,P] = corrcoef(Growth_values_beta,Loss_Indifference); % No

% Is there an association with Interaction Factor in Earnings with IDs?

Interaction_values_beta = betas_first(4,:)'; % Betas associated with Interaction factor

[R,P] = corrcoef(Interaction_values_beta,SevenDown); % No
[R,P] = corrcoef(Interaction_values_beta,SevenUp); % No
[R,P] = corrcoef(Interaction_values_beta,AADIS_data_final); % Nein
[R,P] = corrcoef(Interaction_values_beta,Uncertainty_data_final); % Nyet
[R,P] = corrcoef(Interaction_values_beta,save_risk); % No
[R,P] = corrcoef(Interaction_values_beta,save_loss); %No

figure
scatter(save_loss,Interaction_values_beta, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Loss Aversion', 'FontSize', 16);
ylabel  ('Interaction Beta-Earnings Model', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'StayGo_Exploratory_Loss.svg')


%% Logistic regression


% DataUse explain

%DataUse(:,1) % EQ,
%DataUse(:,2) % PNR,
%DataUse(:,3) % DG Fair
%DataUse(:,4) % UG Fair
%DataUse(:,5) % UG Ask
%DataUse(:,6) % DG Ask
%DataUse(:,7) % UG Optimal
%DataUse(:,8) % DG Optimal
%DataUse(:,9) % Mach
%

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
    
    Growth_Total = [Growth_Total; Growth];
    Final_Value_Total = [Final_Value_Total; Final_Value];
    Turn_Left_Total = [Turn_Left_Total; Turn_Left];
    Earnings_Total = [Earnings_Total; Earnings];
    
end

% [B,DEV,STATS] = mnrfit(DataUse(:,1), DataUse(:,6))
% [B,DEV,STATS] = mnrfit(DataUse(:,9), DataUse(:,5))
First_Turn = [];
Second_Turn = [];
Third_Turn = [];
Fourth_Turn = [];
Fifth_Turn = [];
Sixth_Turn = [];
Seventh_Turn = [];
Eighth_Turn = [];
Ninth_Turn = [];
Tenth_Turn = [];

for ii = 1:length(Turn_Left_Total)
    row = round(Turn_Left_Total(ii));
    if row == 1
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        First_Turn = [First_Turn; row];
        
    end
    if row == 2
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Second_Turn = [Second_Turn; row];
        
    end
    
    if row == 3
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Third_Turn = [Third_Turn; row];
        
    end
    
    if row == 4
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Fourth_Turn = [Fourth_Turn; row];
        
    end
    
    if row == 5
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Fifth_Turn = [Fifth_Turn; row];
        
    end
    
    if row == 6
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Sixth_Turn = [Sixth_Turn; row];
        
    end
    
    if row == 7
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Seventh_Turn = [Seventh_Turn; row];
        
    end
    
    if row == 8
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Eighth_Turn = [Eighth_Turn; row];
        
    end
    
    if row == 9
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Ninth_Turn = [Ninth_Turn; row];
        
    end
    
    if row == 10
        row = [Growth_Total(ii), Final_Value_Total(ii), Turn_Left_Total(ii), Earnings_Total(ii)];
        Tenth_Turn = [Tenth_Turn; row];
        
    end
end

%% Okay, we need to rethink this from a different perspective.

% Let's start with the first turn. We have those (0) who left on the first
% turn. Then those who (1) stayed through the second turn.


%% Second turn

% Let's start with people who left on the second turn.

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2)];
% x1 = matrix(:,1);
% x2 = matrix(:,2);
% x = [ones(size(x1)) x1 x2 x1.*x2];
First_Turn(:,3) = 0; %
Second_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3)];
sp = ordinal(y);

[Betas_Second,DEV,STATS_Second] = mnrfit(x,sp);


%% Third Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3)];
sp = categorical(y);

[B_First,DEV,STATS_First] = mnrfit(x,sp);

[B,DEV,STATS] = mnrfit(x,sp);

sp = ordinal(y);

[Betas_Third,DEV,STATS_Third] = mnrfit(x,sp);

%% Fourth Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2); Fourth_Turn(:,1), Fourth_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 0;
Fourth_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3);Fourth_Turn(:,3)];
sp = ordinal(y);
[Betas_Fourth,DEV,STATS_Fourth] = mnrfit(x,sp);

%% Fifth Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2); Fourth_Turn(:,1), Fourth_Turn(:,2); Fifth_Turn(:,1), Fifth_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 0;
Fourth_Turn (:,3) = 0;
Fifth_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3);Fourth_Turn(:,3);Fifth_Turn(:,3)];
sp = ordinal(y);
[Betas_Fifth,DEV,STATS_Fifth] = mnrfit(x,sp);

%% Sixth Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2); Fourth_Turn(:,1), Fourth_Turn(:,2); Fifth_Turn(:,1), Fifth_Turn(:,2); Sixth_Turn(:,1), Sixth_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 0;
Fourth_Turn (:,3) = 0;
Fifth_Turn (:,3) = 0;
Sixth_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3);Fourth_Turn(:,3);Fifth_Turn(:,3);Sixth_Turn(:,3)];
sp = ordinal(y);
[Betas_Sixth,DEV,STATS_Sixth] = mnrfit(x,sp);

%% Seventh Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2); Fourth_Turn(:,1), Fourth_Turn(:,2); Fifth_Turn(:,1), Fifth_Turn(:,2); Sixth_Turn(:,1), Sixth_Turn(:,2); Seventh_Turn(:,1), Seventh_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 0;
Fourth_Turn (:,3) = 0;
Fifth_Turn (:,3) = 0;
Sixth_Turn (:,3) = 0;
Seventh_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3);Fourth_Turn(:,3);Fifth_Turn(:,3);Sixth_Turn(:,3);Seventh_Turn(:,3)];
sp = ordinal(y);
[Betas_Seventh,DEV,STATS_Seventh] = mnrfit(x,sp);

%% Eighth Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2); Fourth_Turn(:,1), Fourth_Turn(:,2); Fifth_Turn(:,1), Fifth_Turn(:,2); Sixth_Turn(:,1), Sixth_Turn(:,2); Seventh_Turn(:,1), Seventh_Turn(:,2); Eighth_Turn(:,1), Eighth_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 0;
Fourth_Turn (:,3) = 0;
Fifth_Turn (:,3) = 0;
Sixth_Turn (:,3) = 0;
Seventh_Turn (:,3) = 0;
Eighth_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3);Fourth_Turn(:,3);Fifth_Turn(:,3);Sixth_Turn(:,3);Seventh_Turn(:,3);Eighth_Turn(:,3)];
sp = ordinal(y);
[Betas_Eighth,DEV,STATS_Eighth] = mnrfit(x,sp);

%% Ninth Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2); Fourth_Turn(:,1), Fourth_Turn(:,2); Fifth_Turn(:,1), Fifth_Turn(:,2); Sixth_Turn(:,1), Sixth_Turn(:,2); Seventh_Turn(:,1), Seventh_Turn(:,2); Eighth_Turn(:,1), Eighth_Turn(:,2); Ninth_Turn(:,1), Ninth_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 0;
Fourth_Turn (:,3) = 0;
Fifth_Turn (:,3) = 0;
Sixth_Turn (:,3) = 0;
Seventh_Turn (:,3) = 0;
Eighth_Turn (:,3) = 0;
Ninth_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3);Fourth_Turn(:,3);Fifth_Turn(:,3);Sixth_Turn(:,3);Seventh_Turn(:,3);Eighth_Turn(:,3);Ninth_Turn(:,3)];
sp = ordinal(y);
[Betas_Ninth,DEV,STATS_Ninth] = mnrfit(x,sp);

%% Tenth Turn?

x = [First_Turn(:,1), First_Turn(:,2); Second_Turn(:,1), Second_Turn(:,2); Third_Turn(:,1), Third_Turn(:,2); Fourth_Turn(:,1), Fourth_Turn(:,2); Fifth_Turn(:,1), Fifth_Turn(:,2); Sixth_Turn(:,1), Sixth_Turn(:,2); Seventh_Turn(:,1), Seventh_Turn(:,2); Eighth_Turn(:,1), Eighth_Turn(:,2); Ninth_Turn(:,1), Ninth_Turn(:,2); Tenth_Turn(:,1), Tenth_Turn(:,2)];
First_Turn(:,3) = 0;
Second_Turn (:,3) = 0;
Third_Turn (:,3) = 0;
Fourth_Turn (:,3) = 0;
Fifth_Turn (:,3) = 0;
Sixth_Turn (:,3) = 0;
Seventh_Turn (:,3) = 0;
Eighth_Turn (:,3) = 0;
Ninth_Turn (:,3) = 0;
Tenth_Turn (:,3) = 1;
y = [First_Turn(:,3);Second_Turn(:,3);Third_Turn(:,3);Fourth_Turn(:,3);Fifth_Turn(:,3);Sixth_Turn(:,3);Seventh_Turn(:,3);Eighth_Turn(:,3);Ninth_Turn(:,3);Tenth_Turn(:,3)];
sp = ordinal(y);
[Betas_Tenth,DEV,STATS_Tenth] = mnrfit(x,sp);

%% Bar plot Growth


data = [Betas_Second, Betas_Third, Betas_Fourth, Betas_Fifth, Betas_Sixth, Betas_Seventh, Betas_Eighth, Betas_Ninth, Betas_Tenth];
x = [linspace(1,9,9); linspace(1,9,9); linspace(1,9,9)];
B1Er = [STATS_Second.se(1)*2,STATS_Third.se(1)*2,STATS_Fourth.se(1)*2,STATS_Fifth.se(1)*2,STATS_Sixth.se(1)*2,STATS_Seventh.se(1)*2,STATS_Eighth.se(1)*2,STATS_Ninth.se(1)*2,STATS_Tenth.se(1)*2];
B2Er = [STATS_Second.se(2)*2,STATS_Third.se(2)*2,STATS_Fourth.se(2)*2,STATS_Fifth.se(2)*2,STATS_Sixth.se(2)*2,STATS_Seventh.se(2)*2,STATS_Eighth.se(2)*2,STATS_Ninth.se(2)*2,STATS_Tenth.se(2)*2];
B3Er = [STATS_Second.se(3)*2,STATS_Third.se(3)*2,STATS_Fourth.se(3)*2,STATS_Fifth.se(3)*2,STATS_Sixth.se(3)*2,STATS_Seventh.se(3)*2,STATS_Eighth.se(3)*2,STATS_Ninth.se(3)*2,STATS_Tenth.se(3)*2];
err = [B1Er;B2Er;B3Er];

figure
bar(x(1,:),data(1,:))
ax = gca;
ax.FontSize = 12;
xlim ([0 10]);
box off
xlabel ('Betas', 'FontSize', 16);
ylabel  ('Beta Weight', 'FontSize', 16);
set(gcf,'color','w');
set(gca, 'XTick', 1:9, 'XTickLabels', {'Turn 2 ','Turn 3 ','Turn 4 ','Turn 5 ','Turn 6 ','Turn 7 ','Turn 8 ','Turn 9 ','Turn 10 '})
title('Logistic regressions of Growth Factors')

hold on

% Standard Error

er = errorbar(x(1,:),data(1,:),err(1,:)); %errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
hold off


saveas(gcf,'Logistic_regression_Growth.svg')

%% Bar plot Final Values


data = [Betas_Second, Betas_Third, Betas_Fourth, Betas_Fifth, Betas_Sixth, Betas_Seventh, Betas_Eighth, Betas_Ninth, Betas_Tenth];
x = [linspace(1,9,9); linspace(1,9,9); linspace(1,9,9)];
B1Er = [STATS_Second.se(1)*2,STATS_Third.se(1)*2,STATS_Fourth.se(1)*2,STATS_Fifth.se(1)*2,STATS_Sixth.se(1)*2,STATS_Seventh.se(1)*2,STATS_Eighth.se(1)*2,STATS_Ninth.se(1)*2,STATS_Tenth.se(1)*2];
B2Er = [STATS_Second.se(2)*2,STATS_Third.se(2)*2,STATS_Fourth.se(2)*2,STATS_Fifth.se(2)*2,STATS_Sixth.se(2)*2,STATS_Seventh.se(2)*2,STATS_Eighth.se(2)*2,STATS_Ninth.se(2)*2,STATS_Tenth.se(2)*2];
B3Er = [STATS_Second.se(3)*2,STATS_Third.se(3)*2,STATS_Fourth.se(3)*2,STATS_Fifth.se(3)*2,STATS_Sixth.se(3)*2,STATS_Seventh.se(3)*2,STATS_Eighth.se(3)*2,STATS_Ninth.se(3)*2,STATS_Tenth.se(3)*2];
err = [B1Er;B2Er;B3Er];

figure
bar(x(2,:),data(2,:))
ax = gca;
ax.FontSize = 12;
xlim ([0 10]);
box off
xlabel ('Betas', 'FontSize', 16);
ylabel  ('Beta Weight', 'FontSize', 16);
set(gcf,'color','w');
set(gca, 'XTick', 1:9, 'XTickLabels', {'Turn 2 ','Turn 3 ','Turn 4 ','Turn 5 ','Turn 6 ','Turn 7 ','Turn 8 ','Turn 9 ','Turn 10 '})
title('Logistic regressions of Final Values Factors')

hold on

% Standard Error

er = errorbar(x(2,:),data(2,:),err(2,:)); %errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
hold off


saveas(gcf,'Logistic_regression_FinalValues.svg')

%% Bar plot Interaction


data = [Betas_Second, Betas_Third, Betas_Fourth, Betas_Fifth, Betas_Sixth, Betas_Seventh, Betas_Eighth, Betas_Ninth, Betas_Tenth];
x = [linspace(1,9,9); linspace(1,9,9); linspace(1,9,9)];
B1Er = [STATS_Second.se(1)*2,STATS_Third.se(1)*2,STATS_Fourth.se(1)*2,STATS_Fifth.se(1)*2,STATS_Sixth.se(1)*2,STATS_Seventh.se(1)*2,STATS_Eighth.se(1)*2,STATS_Ninth.se(1)*2,STATS_Tenth.se(1)*2];
B2Er = [STATS_Second.se(2)*2,STATS_Third.se(2)*2,STATS_Fourth.se(2)*2,STATS_Fifth.se(2)*2,STATS_Sixth.se(2)*2,STATS_Seventh.se(2)*2,STATS_Eighth.se(2)*2,STATS_Ninth.se(2)*2,STATS_Tenth.se(2)*2];
B3Er = [STATS_Second.se(3)*2,STATS_Third.se(3)*2,STATS_Fourth.se(3)*2,STATS_Fifth.se(3)*2,STATS_Sixth.se(3)*2,STATS_Seventh.se(3)*2,STATS_Eighth.se(3)*2,STATS_Ninth.se(3)*2,STATS_Tenth.se(3)*2];
err = [B1Er;B2Er;B3Er];

figure
bar(x(3,:),data(3,:))
ax = gca;
ax.FontSize = 12;
xlim ([0 10]);
box off
xlabel ('Betas', 'FontSize', 16);
ylabel  ('Beta Weight', 'FontSize', 16);
set(gcf,'color','w');
set(gca, 'XTick', 1:9, 'XTickLabels', {'Turn 2 ','Turn 3 ','Turn 4 ','Turn 5 ','Turn 6 ','Turn 7 ','Turn 8 ','Turn 9 ','Turn 10 '})
title('Logistic regressions of Interaction Factors')

hold on

% Standard Error

er = errorbar(x(3,:),data(3,:),err(3,:)); %errorbar(x,data,errlow,errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
hold off


saveas(gcf,'Logistic_regression_Interaction.svg')


