%% Initialization

clear all
close all
clc

% This script analyzes the EQ Thesis data. 

% Daniel Sazhin
% 05/03/19
% Thesis

%% Useful things for the FUTURE

% Histogram
% Bar Plot (with error bars)
% Scatterplots
% Logistic Regression Analysis
% Correlation Analysis

load Behavioral_Data.mat
load Behavioral_Survey_Data.mat 
load Survey_Data.mat 

%% DG_A Earnings 

% Total earnings are accumulated over all trials.

% For DG, it is simply the means * number of trials.

DG_A_Earnings = (1-DG_A_Mean_Offers)*25;

%% UG_B Acceptance Earnings

load Behavioral_Data.mat

X1 = cell2mat(X1(2:end,:)); % Convert the data into a matrix file.
[m,n]= size(X1); % Find the width of the data.
X1_Subject_Earnings = [];
 
for ii = 1:length(X1) % For all the subjects
    NumberAccepted = ((sum(X1(ii,:)))- n); % Find out the number of offers accepted. 
    saveme = NumberAccepted*.1; % Multiply by the amount earned.
    X1_Subject_Earnings = [X1_Subject_Earnings, saveme]; % Save it. 
end    

X2 = cell2mat(X2(2:end,:));
[m,n]= size(X2);
X2_Subject_Earnings = [];

for ii = 1:length(X2)
    NumberAccepted = ((sum(X2(ii,:)))- n);
    saveme = NumberAccepted*.2;
    X2_Subject_Earnings = [X2_Subject_Earnings, saveme];
end    

X3 = cell2mat(X3(2:end,:));
[m,n]= size(X3);
X3_Subject_Earnings = [];

for ii = 1:length(X3)
    NumberAccepted = ((sum(X3(ii,:)))- n);
    saveme = NumberAccepted*.3;
    X3_Subject_Earnings = [X3_Subject_Earnings, saveme];
end    

X4 = cell2mat(X4(2:end,:));
[m,n]= size(X4);
X4_Subject_Earnings = [];

for ii = 1:length(X4)
    NumberAccepted = ((sum(X4(ii,:)))- n);
    saveme = NumberAccepted*.4;
    X4_Subject_Earnings = [X4_Subject_Earnings, saveme];
end    

X5 = cell2mat(X5(2:end,:));
[m,n]= size(X5);
X5_Subject_Earnings = [];

for ii = 1:length(X5)
    NumberAccepted = ((sum(X5(ii,:)))- n);
    saveme = NumberAccepted*.5;
    X5_Subject_Earnings = [X5_Subject_Earnings, saveme];
end    

X6 = cell2mat(X6(2:end,:));
[m,n]= size(X6);
X6_Subject_Earnings = [];

for ii = 1:length(X6)
    NumberAccepted = ((sum(X6(ii,:)))- n);
    saveme = NumberAccepted*.6;
    X6_Subject_Earnings = [X6_Subject_Earnings, saveme];
end   

UG_B_Earnings = X1_Subject_Earnings + X2_Subject_Earnings + X3_Subject_Earnings + X4_Subject_Earnings + X5_Subject_Earnings + X6_Subject_Earnings;

%% UG_A Earnings

% This is less straightforward. Assume refusals below a given threshold.

ExpectedEarnings = 1- UG_A;
UG_A_ExpectedEarnings = sum(ExpectedEarnings');

% hist(UG_A(:) ./ length(UG_A))
%hist(UG_A(:))./ length(UG_A)
% Find proportion of acceptance/rejections.

load Behavioral_Data.mat

X0 = cell2mat(X0(2:end,:)); % Convert the data into a matrix file.
X1 = cell2mat(X1(2:end,:));
X2 = cell2mat(X2(2:end,:));
X3 = cell2mat(X3(2:end,:));
X4 = cell2mat(X4(2:end,:));
X5 = cell2mat(X5(2:end,:));
X6 = cell2mat(X6(2:end,:));

[TotalSubjects,n] = size(X0);

% Find proportion of acceptance/rejections in our subjects

Reject = (sum(X0(:) == 1)); % Number of rejections
Accept = (sum(X0(:) == 2));
Total = Accept + Reject;
X0_Accept = (Accept/Total);
VarX0 = n * X0_Accept * (1 - X0_Accept)

Reject = (sum(X1(:) == 1)); % Number of rejections
Accept = (sum(X1(:) == 2));
Total = Accept + Reject;
X1_Accept = (Accept/Total);
VarX1 = n * X1_Accept * (1 - X1_Accept)

Reject = (sum(X2(:) == 1)); % Number of rejections
Accept = (sum(X2(:) == 2));
Total = Accept + Reject;
X2_Accept = (Accept/Total);
VarX2 = n * X2_Accept * (1 - X2_Accept)

Reject = (sum(X3(:) == 1)); % Number of rejections
Accept = (sum(X3(:) == 2));
Total = Accept + Reject;
X3_Accept = (Accept/Total);
VarX3 = n * X3_Accept * (1 - X3_Accept)

Reject = (sum(X4(:) == 1)); % Number of rejections
Accept = (sum(X4(:) == 2));
Total = Accept + Reject;
X4_Accept = (Accept/Total);
VarX4 = n * X4_Accept * (1 - X4_Accept)

Reject = (sum(X5(:) == 1)); % Number of rejections
Accept = (sum(X5(:) == 2));
Total = Accept + Reject;
X5_Accept = (Accept/Total);
VarX5 = n * X5_Accept * (1 - X5_Accept)

Reject = (sum(X6(:) == 1)); % Number of rejections
Accept = (sum(X6(:) == 2));
Total = Accept + Reject;
X6_Accept = (Accept/Total);
VarX6 = n * X6_Accept * (1 - X6_Accept)

UG_A_Earnings = [];

for ii = 1:length(UG_A)
    subject = UG_A(ii,:);
    
    % Count up each offer they have made.
    % Multiply it by the acceptance rate.
    % Multiply by offer amount.
    
    X0_offer = (sum(subject(:) == 0))*X0_Accept*1;
    X1_offer = (sum(subject(:) == .1))*X1_Accept*.9;
    X2_offer = (sum(subject(:) == .2))*X2_Accept*.8;
    X3_offer = (sum(subject(:) == .3))*X3_Accept*.7;
    X4_offer = (sum(subject(:) == .4))*X4_Accept*.6;
    X5_offer = (sum(subject(:) == .5))*X5_Accept*.5;
    X6_offer = (sum(subject(:) == .6))*X6_Accept*.4;
    X7_offer = (sum(subject(:) == .7))*.3; % Assume always accept
    X8_offer = (sum(subject(:) == .8))*.2;
    X9_offer = (sum(subject(:) == .9))*.1;
    X10_offer = (sum(subject(:) == 1))*0;
    
    SubjectEarnings = X0_offer+X1_offer+X2_offer+X3_offer+X4_offer+X5_offer+X6_offer+X7_offer+X8_offer+X9_offer+X10_offer;
    UG_A_Earnings = [UG_A_Earnings,SubjectEarnings];
end    
    

%% Total Hypothetical Earnings

TotalEarnings = UG_A_Earnings + UG_B_Earnings + DG_A_Earnings;

%% Hypothesis Testing
 
% Hypothesis

% Higher EQ leads to higher earnings due to less rejection of offers as
% Player A

[R,P] = corrcoef(UG_A_Earnings,TotalEQScore);

subplot(2,4,1)
scatter(UG_A_Earnings, TotalEQScore)
title('UG_A Earnings and Total EQ Nsig')
% False

% Hypothesis:

% Higher EQ leads to greater rejection of unfair offers in UG.
% Unfair offers are classifed as $0, $1, $2, $3.

% Proportion of rejection of unfair offers.

% Group together all unfair offers.

load Behavioral_Data.mat

Unfair_Offers = [X0,X1,X2,X3];
Unfair_Offers = cell2mat(Unfair_Offers(2:end,:)); % Convert the data into a matrix file.

Proportion = [];
for ii = 1:length(Unfair_Offers)
    subject = Unfair_Offers(ii,:);
    Reject = (sum(subject(:) == 1));
    Accept = (sum(subject(:) == 2));
    saveme = Reject/(Accept+Reject);
    Proportion = [Proportion,saveme];
end    

% Total_Offers = [X0,X1,X2,X3,X4,X5,X6];
% Total_Offers = cell2mat(Total_Offers(2:end,:)); % Convert the data into a matrix file.
% 
% Proportion = [];
% for ii = 1:length(Total_Offers)
%     subject = Total_Offers(ii,:);
%     Reject = (sum(subject(:) == 1));
%     Accept = (sum(subject(:) == 2));
%     saveme = Reject/(Accept+Reject);
%     Proportion = [Proportion,saveme];
% end   

[R,P] = corrcoef(Proportion,TotalEQScore);

subplot(2,4,2)
scatter(Proportion, TotalEQScore)
title 'Proportion and EQ Nsig'

% Hypothesis 

% Greater Negative PNR leads to greater rejections.

[R,P] = corrcoef(Proportion,PNRScore);

subplot(2,4,3)
scatter(Proportion, PNRScore)
title 'Proportion and PNR Sig'

% Significant Result.

% Hypothesis 

% Higher Mach and Higher Self Control lead to higher total earnings.

[R,P] = corrcoef(PNRScore, Overall_Mach_Score)
subplot(2,4,4)
scatter(PNRScore, Overall_Mach_Score)
title 'PNRScore, Overall_Mach_Score'

% Hypothesis 
% Higher Emotionality leads to more rejections of unfair offers.

[R,P] = corrcoef(EmotionalityScore, Proportion)
subplot(2,4,5)
scatter(EmotionalityScore, Proportion)
title 'Proportion and Emotionality NSig'
% false

% Hypothesis 

% Higher Self Control lead to higher total earnings.

[R,P] = corrcoef(TotalEarnings, SelfcontrolScore)
subplot(2,4,6)
scatter(TotalEarnings, SelfcontrolScore)
title 'Total Earnings and Self Control Sig'
% Significant

% Hypothesis 

% Higher EQ lead to higher total earnings.

[R,P] = corrcoef(TotalEarnings, TotalEQScore)
subplot(2,4,7)
scatter(TotalEarnings, TotalEQScore)
title 'Total Earnings and Total EQ Sig')

% Significant

% Why is this?

% Because higher DG_A Offers:

[R,P] = corrcoef(DG_A_Earnings, TotalEQScore)
subplot(2,4,8)
scatter(DG_A_Earnings, TotalEQScore)
title 'DG_A, Total EQ Sig'

figure
hist(UG_A(:))
title 'Histogram of Ultimatum Game Offers as Player A'
xlabel 'Offers'
ylabel 'Frequency'

figure
hist(DG_A(:))
title 'Histogram of Dictator Game Offers as Player A'
xlabel 'Offers'
ylabel 'Frequency'

%% UG_B

% What is the proportion of rejections for each class of offer?

figure
plot([X0_Accept,X1_Accept,X2_Accept,X3_Accept,X4_Accept,X5_Accept,X6_Accept])

%% Fair/Selfish Test

Fair = cell2mat(DG_Ask(2:end,:));

test = [Fair,TotalEQScore'];
Selfish = [];
Fair2= [];

for ii = 1:length(Fair)
    row = test(ii,:);
    if row(1) == 1
        saveme = test(ii,:);
        Fair2 = [Fair2;saveme];
    end
    if row(1) == 2
        saveas2 = test(ii,:);
        Selfish = [Selfish;saveas2];
    end
end    
    
[H,P,CI,STATS] = ttest2(Fair2(:,2), Selfish(:,2));

% No signficant difference between samples.

% No difference in EQ between self-reported Selfish and Fair participants.)

%% Multiple regression 

% Higher Mach and Higher Self Control lead to higher total earnings.

% Multiple regression interaction matlab. 

x1 = Overall_Mach_Score';
x2 = TotalEQScore';    % Contains NaN data
y = TotalEarnings';

X = [ones(size(x1)) x1 x2 x1.*x2];
b = regress(y,X)    % Removes NaN data

[~,~,r,rint] = regress(y,X,0.05)
contain0 = (rint(:,1)<0 & rint(:,2)>0);
idx = find(contain0==false)
figure
hold on
scatter(y,r)
scatter(y(idx),r(idx),'b','filled')
xlabel("Mach and EQ")
ylabel("Earnings")
hold off

%% Correlations

% R = corrcoef(A,B)

[R,P] = corrcoef(TotalEQScore,TotalEarnings);

% High Mach, Highself Control, better earnings

% %%
% 
% % clearvars -except UG_A_Mean_Offers DG_A_Mean_Offers
% 
% [R,P,RL,RU] = corrcoef(UG_A_Mean_Offers, TotalEQScore);
% 
% data1 = UG_A_Mean_Offers;
% data2 =  TotalEQScore;
% 
% figure; 
% scatter(data1, data2) %, '+', 'MarkerFaceColor', 'k');
% ylabel('TotalEQScore');
% xlabel ('UG_A_Mean_Offers');
% box 'on'
% % ticks off
% set(gca,'Ticklength',[0 0])
% axis square;
% r = corrcoef(data1, data2);
% disp(r(1,2));
% str=['r= ',num2str(r(1,2))]
% T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
% set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%% Correlations between surveys

[R,P] = corrcoef(PNRScore,TotalEQScore)

figure
subplot(2,4,1)
scatter(PNRScore,TotalEQScore)
title 'PNRScore,TotalEQScore'

[R,P] = corrcoef(PNRScore,Overall_Mach_Score)

subplot(2,4,2)
scatter(PNRScore,Overall_Mach_Score)
title 'PNRScore,Overall_Mach_Score'

[R,P] = corrcoef(TotalEQScore,Overall_Mach_Score)

subplot(2,4,3)
scatter(TotalEQScore,Overall_Mach_Score)
title 'TotalEQScore,Overall_Mach_Score'

[R,P] = corrcoef(Fair,Overall_Mach_Score)
subplot(2,4,4)
scatter(Proportion,Overall_Mach_Score)
title 'Proportion,Overall_Mach_Score'

[R,P] = corrcoef(UG_B_Earnings,PNRScore)
subplot(2,4,5)
scatter(UG_B_Earnings,PNRScore)
title 'UG_B_Earnings,PNRScore'

%% T-tests

meanoffersA=(cell2mat(DG_Fair(2:end,:)));
meanoffersB=(cell2mat(UG_Fair(2:end,:)));
 meanoffersC=(cell2mat(UG_Ask(2:end,:)));
 sum(meanoffersC(:) == 1); % 1 = fair;
 meanoffersD=(cell2mat(DG_Ask(2:end,:)));
  [H,P,CI,STATS] = ttest2(mean(DG_A),mean(UG_A));
 sum(meanoffersD(:) == 1) ;
 
 meanoffersE=(cell2mat(UG_Optimal(2:end,:)));
 meanoffersF=(cell2mat(DG_Optimal(2:end,:)));
 [H,P,CI,STATS] = ttest(meanoffersE,meanoffersF);
 
 %% Psychometric
 
 figure
hist(TotalEQScore)
title('Emotional Intelligence')
xlabel('Level of EI')
ylabel('Frequency')

figure
hist(PNRScore)
title('Personal Norms of Reciprocity')
xlabel('Level of PNR')
ylabel('Frequency')

figure
hist(Overall_Mach_Score)
title('Machiavellianism')
xlabel('Level of Machiavellianism')
ylabel('Frequency')

figure
scatter(Proportion, SelfcontrolScore);
 

%% Test for order effects

Order_test = [UG_A_Earnings',UG_B_Earnings',DG_A_Earnings',Proportion',TotalEQScore'];

% Sort Order_test by groups.

% Find subjects by row.

UG_First_Subjects = [];
DG_First_Subjects = [];

for ii = 1:length(UG_First)
    subject = UG_First(ii);
    saveme = Order_test(subject,:);
    UG_First_Subjects = [UG_First_Subjects; saveme];
end    

for ii = 1:length(DG_First)
    subject = DG_First(ii);
    saveme = Order_test(subject,:);
    DG_First_Subjects = [DG_First_Subjects; saveme];
end   

% Test if earnings, EQ are different between groups.

[H,P_UG_A_Earnings,CI,STATS] = ttest2(UG_First_Subjects(:,1), DG_First_Subjects(:,1));
[H,P_UG_B_Earnings,CI,STATS] = ttest2(UG_First_Subjects(:,2), DG_First_Subjects(:,2));
[H,P_DG_A_Earnings,CI,STATS] = ttest2(UG_First_Subjects(:,3), DG_First_Subjects(:,3));
[H,P_Proportion_Rejections,CI,STATS] = ttest2(UG_First_Subjects(:,4), DG_First_Subjects(:,4));
[H,P_EQ_Rejections,CI,STATS] = ttest2(UG_First_Subjects(:,5), DG_First_Subjects(:,5));

[R,P] = corrcoef(UG_First_Subjects(:,4),UG_First_Subjects(:,5)); % Correlation of EQ and Proportion reject
; % Correlation of EQ and DG_A_Earnings


% What about proportion of max to min? Scaled if you will.

% Max in UG_A 13.2903
% Max in DG_A 25.00

% How much did each subject make relative to this max?


% Those who play DG_A well, don't tend to play UG_A well. 

Proportion_UGA_Earnings = UG_A_Earnings ./ max(UG_A_Earnings); % 13.2903 max
Proportion_DGA_Earnings = DG_A_Earnings ./ max(DG_A_Earnings); % 25 max
[R,P] = corrcoef(Proportion_UGA_Earnings,Proportion_DGA_Earnings)

figure
scatter(Proportion_UGA_Earnings, Proportion_DGA_Earnings)
xlabel('UGA')
ylabel('DGA')
axis([0 1 0 1])
title('Proportion of Earnings')

Both = (Proportion_UGA_Earnings + Proportion_DGA_Earnings) ./ 2;

[R,P] = corrcoef(Both, TotalEQScore)

figure
scatter(Both, TotalEQScore)
xlabel('Both')
ylabel('TotalEQScore')
title('Proportion of Earnings for Both and EQ Score')

% Those who are the most strategic, have the most EQ.


% What is the EI of the generous folks in the DG Condition?

% Define generous as keeping less than 17.5 (.3 * 25)

Generous = [DG_A_Earnings',TotalEQScore',PNRScore',Overall_Mach_Score'];
GenerousSubjects = [];
SelfishSubjects = [];
for ii = 1:length(Generous)
    row = Generous(ii,:);
    if row(1,1) <= 17.5
       GenerousSubjects = [GenerousSubjects;row];
    else
        SelfishSubjects = [SelfishSubjects;row];
    end
end    

[R,P] = corrcoef(GenerousSubjects(:,1), GenerousSubjects(:,2))

figure
scatter(GenerousSubjects(:,1), GenerousSubjects(:,2))
xlabel('How Selfish Were Subjects')
ylabel('EI')

title('Generosity and EI')

% EI does not predict

[~,P] = corrcoef(SelfishSubjects(:,1), SelfishSubjects(:,2))

% Strength of .5 offer irrespective of EI. How engrained this is!


% Omit .5 offers in UG. What happens?

UGLess = [mean(UG_A')',TotalEQScore',PNRScore',Overall_Mach_Score'];
Lessthanhalf =  [];
Morethanhalf = [];
for ii = 1:length(UGLess)
    row = UGLess(ii,:);
    if row(1,1) < .5
       Lessthanhalf = [Lessthanhalf;row];
    else
        Morethanhalf = [Morethanhalf;row];
    end
end  

[R,P] = corrcoef(Lessthanhalf(:,1), Lessthanhalf(:,2))

figure
scatter(Morethanhalf(:,1), Morethanhalf(:,2))
xlabel('UG Offer')
ylabel('EI')

title('Generosity and EI')

% DG/UG consistent shows people are pretty consistent. A lot in UG, offer a
% lot in DG. 

% possible limitation. People thought games are the same. Due to
% correlation. 


%% [R,P] = corrcoef(PNRScore, Overall_Mach_Score)

figure
scatter(UG_B_Earnings, TotalEQScore)
title 'UG_B_Earnings, Overall_Mach_Score'
[R,P] = corrcoef(UG_B_Earnings,TotalEQScore)

figure
scatter(UG_A_Earnings,TotalEQScore)
title 'UG_A_Earnings, EI Score'


MedianSplit = [TotalEQScore', TotalEarnings', meanoffersA, meanoffersB, meanoffersC, meanoffersD, meanoffersE, meanoffersF, Overall_Mach_Score'];
HighEI = [];
LowEI = [];
for ii = 1:length(MedianSplit)
    row = MedianSplit(ii,:);
    if row(1,9) > 71;
       HighEI = [HighEI; MedianSplit(ii,:)];
    end   
    if row(1,9) <= 71
       LowEI = [LowEI; MedianSplit(ii,:)];
    end   
end   

HighEI(:,3) % DG Fair A 
HighEI(:,4) % UG Fair B
HighEI(:,5) % UG Ask C
HighEI(:,6) % DG Ask D 
HighEI(:,7) % UG Optimal E
HighEI(:,8) % DG Optimal F

figure
scatter(meanoffersF,TotalEQScore)

%%

LowEIselfish = [];
for ii = 1:length(LowEI)
    row = LowEI(ii,6)
        if row == 1
            LowEIselfish = [LowEIselfish, row];
        end
end

HighEIselfish = [];
for ii = 1:length(HighEI)
    row = HighEI(ii,6)
        if row == 1
            HighEIselfish = [HighEIselfish, row];
        end
end

% 53 High vs 56 self report Fair in UG

% 30 High report fair in DG 46 report fair as low EI

%meanoffersC=(cell2mat(UG_Ask(2:end,:)));
% sum(meanoffersC(:) == 1); % 1 = fair;

%% Logistic regression

DataUse = MedianSplit;

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

[B,DEV,STATS] = mnrfit(DataUse(:,1), DataUse(:,6))
[B,DEV,STATS] = mnrfit(DataUse(:,9), DataUse(:,5))

%% Scatterplots prototype

% [R,P] = corrcoef(TotalEarnings, TotalEQScore)
% figure
% scatter(TotalEarnings, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% ax = gca
% ax.FontSize = 12
% xlabel ('Total Earnings in Dollars', 'FontSize', 16);
% ylabel  ('EI Score', 'FontSize', 16);
%i = lsline;
%i.LineWidth = 5;
%i.Color = [0 0 0];
% set(gcf,'color','w');

%% Histogram prototype

% figure
% h = histogram(DG_A(:));
% counts = h.Values;
% h.NumBins = 11
% ax = gca
% ax.FontSize = 12
% xlabel ('Offer Amount','FontSize', 16)
% ylabel ('Frequency','FontSize', 16)
% set(gca,'box','off')
% set(gcf,'color','w');

%% Figure 1

figure
h = histogram(TotalEQScore);
counts = h.Values;
h.NumBins = 12
ax = gca
ax.FontSize = 12
xlabel ('Emotional Intelligence','FontSize', 16)
ylabel ('Frequency','FontSize', 16)
set(gca,'box','off')
set(gcf,'color','w');

saveas(gcf,'Fig1.png')

%% Figure 3

figure
h = histogram(UG_A(:));
counts = h.Values;
h.NumBins = 11
ax = gca
ax.FontSize = 9
xlabel ('Ultimatum Game Offers as Player A','FontSize', 16)
ylabel ('Frequency','FontSize', 16)
set(gca,'box','off')
set(gcf,'color','w');

saveas(gcf,'Fig3.png')

%% Figure 6

figure
h = histogram(DG_A(:));
counts = h.Values;
h.NumBins = 11
ax = gca
ax.FontSize = 9
xlabel ('Dictator Game Offers as Player A','FontSize', 16)
ylabel ('Frequency','FontSize', 16)
set(gca,'box','off')
set(gcf,'color','w');

saveas(gcf,'Fig6.png')

%% Scatterplot

[R,P] = corrcoef(TotalEarnings, TotalEQScore)
figure
scatter(TotalEarnings, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Total Earnings in Dollars', 'FontSize', 16);
ylabel  ('EI Score', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'TotalEarnings.png')






%% Figure 2

[R,P] = corrcoef(Proportion, TotalEQScore)
figure
scatter(Proportion, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('P(Rejections)', 'FontSize', 16);
ylabel  ('EI Score', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'Fig2.png')

%% Figure 5

[R,P] = corrcoef(UG_A_Earnings, TotalEQScore);
figure
scatter(UG_A_Earnings, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca;
ax.FontSize = 12;
xlabel ('Total Earnings as Player A in UG Task', 'FontSize', 16);
ylabel  ('EI Score', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'Fig5.png')

%% Figure 7

[R,P] = corrcoef(DG_A_Earnings, TotalEQScore);
figure
scatter(DG_A_Earnings, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5);
ax = gca;
ax.FontSize = 12;
xlabel ('Total Earnings as Player A in DG Task', 'FontSize', 16);
ylabel  ('EI Score', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'Fig7.png')

%% Figure 8

[R,P] = corrcoef(TotalEarnings, TotalEQScore)
figure
scatter(TotalEarnings, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Total Earnings', 'FontSize', 16);
ylabel  ('EI Score', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'Fig8.png')

%% Figure 9

[R,P] = corrcoef(PNRScore, TotalEQScore)
figure
scatter(PNRScore, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('PNR Score', 'FontSize', 16);
ylabel  ('EI Score', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'Fig9.png')

%% Figure 10

[R,P] = corrcoef(Overall_Mach_Score, TotalEQScore)
figure
scatter(Overall_Mach_Score, TotalEQScore, 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('MACH Score', 'FontSize', 16);
ylabel  ('EI Score', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'Fig10.png')

%% Bar Plot

load Behavioral_Data.mat

X0 = cell2mat(X0(2:end,:)); % Convert the data into a matrix file.
X1 = cell2mat(X1(2:end,:));
X2 = cell2mat(X2(2:end,:));
X3 = cell2mat(X3(2:end,:));
X4 = cell2mat(X4(2:end,:));
X5 = cell2mat(X5(2:end,:));
X6 = cell2mat(X6(2:end,:));

X0 = X0 - 1;
X1 = X1 - 1;
X2 = reshape(X2,[],1) - 1;
X3 = reshape(X3,[],1) - 1;
X4 = reshape(X4,[],1) - 1;
X5 = reshape(X5,[],1) - 1;
X6 = reshape(X6,[],1) - 1;

% Standard Error 

X0Er = std(X0) / sqrt(length(X0));
X1Er = std(X1) / sqrt(length(X1));
X2Er = std(X2) / sqrt(length(X2));
X3Er = std(X3) / sqrt(length(X3));
X4Er = std(X4) / sqrt(length(X4));
X5Er = std(X5) / sqrt(length(X5));
X6Er = std(X6) / sqrt(length(X6));

err = [X0Er,X1Er,X2Er,X3Er,X4Er,X5Er,X6Er] * 2;

% What is the proportion of rejections for each class of offer?

figure
x = linspace(0,.6,7)
y = [X0_Accept,X1_Accept,X2_Accept,X3_Accept,X4_Accept,X5_Accept,X6_Accept];
b = bar(x,y)
data = y
errhigh = [X0_Accept+X0Er X1_Accept+X1Er X2_Accept+X2Er X3_Accept+X3Er X4_Accept+X4Er X5_Accept+X5Er X6_Accept+X6Er];
errlow  = [X0_Accept-X0Er X1_Accept-X1Er X2_Accept-X2Er X3_Accept-X3Er X4_Accept-X4Er X5_Accept-X5Er X6_Accept-X6Er];
figure

bar(x,data)    
ax = gca
ax.FontSize = 12
xlim ([-.05 .65])
box off
xlabel ('Offers', 'FontSize', 16);
ylabel  ('P(Accept)', 'FontSize', 16);
set(gcf,'color','w');

hold on

er = errorbar(x,data,err) %errorbar(x,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 1
hold off

saveas(gcf,'Fig4.png')

%% Test for trim

Trimmed = [TotalEQScore',TotalEarnings'];
TrimmedTest = [];
Discard = [];
dontkeepme = [];
for ii = 1:length(Trimmed)
    row = Trimmed(ii,:);
    if row(1,2) > 25;
        keepme = row;
        TrimmedTest = [TrimmedTest;keepme];
    end
    if row(1,2) < 25;
        dontkeepme = row;
    end    
;
    Discard = [Discard; dontkeepme];
end    


%% Figure 11

[R,P] = corrcoef(TrimmedTest(:,1), TrimmedTest(:,2));
figure
scatter(TrimmedTest(:,1), TrimmedTest(:,2), 'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
ax = gca
ax.FontSize = 12
xlabel ('Trimmed EI', 'FontSize', 16);
ylabel  ('Trimmed Earnings', 'FontSize', 16);
i = lsline;
i.LineWidth = 5;
i.Color = [0 0 0];
set(gcf,'color','w');

saveas(gcf,'Fig11.png')