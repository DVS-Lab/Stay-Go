%% Initialization

clear all;
close all;
clc;

[n,t,rawdata] = xlsread('Thesis_Data_0403_leave.csv'); % This converts the CSV file into a matlab file. 

% Rawdata needs to be processed through the exclusion criteria. Then it
% becomes "data"

%% Briefing

% This script will clean up the data from the EQ Qualtrics survey for
% Daniel's thesis.

% Daniel Sazhin
% 03/31/18

% This script takes the CSV file and cleans the data. At the end of it, you
% will have the following useful information from our participants in
% several files:

% Behavioral_Data, Behavioral_Survey_Data, Survey_Data, Demographic_Data

% The actual variables are listed below:

% UG_A_Mean_Offers - Mean offers over all 25 trials as Player A.
% UG_A - Raw data for UG_A results.
% DG_A_Mean_Offers - Mean offers over all 25 trials as Player A.
% UG_B: UG_B_temp has all the data with headers.

% For UG_B, participants accepted or rejected offers. Where 2 is accept and
% 1 is reject.  These have been subdivided into their offer types: 

% X0, X1, X2, X3, X4, X5, X6

% X0 represents the following offer: 10 for Player A, 0 for Player B.

% Behavioral Surveys ----------------------

%   UG_Fair - Fairest offer?
%   UG_Optimal - Offer which gets you the most $$
%   UG_Ask - What was your strategy?

%   DG_Fair - Fairest offer?
%   DG_Optimal - Offer which gets you the most $$
%   DG_Ask - What was your strategy?

% Surveys ---------------------------------

%   PNRScore

%   Overall_Mach_Score
%   Subscales

%       ViewsScore
%       MoralityScore
%       TacticScore

%   TotalEQScore
%   Subscales

%       WellbeingScore
%       SelfcontrolScore
%       EmotionalityScore
%       SociabilityScore

% Demographics ------------------------------------

%   age 

%   gender -------
%       males 
%       females

%   race --------
%       white 
%       black
%       hispanic 
%       indian 
%       asian 
%       PI 
%       mixed
%       other

% education ------------
%       middle_school 
%       high_school 
%       training
%       associate
%       bachelor
%       master 
%       doctor
%       professional

% Use this information toward further analysis.

%% Headers

use = t(1,:); % Pulls out the headers 
headers = cell2table(use); % This turns the headers into a table format
headers = table2array(headers);

%% Exclusion criteria

% 1) If you do not consent, the survey does not let you continue.
% 2) If you are not engaged (1), you get excluded.
% 3) If you do not finish the survey, you get excluded.
% 4) If you do not understand DG/UG Optimal Question, you get excluded.
% 5) If you answer the engagement questions at the end of the survey
% incorrectly, you get excluded.

% Progress

temp = find(strcmp('Finished',rawdata(1,:))); % This finds the Finished column.
Finished = rawdata(:,temp); % Takes out the column from the data table.
Finished = cell2mat(Finished(4:end,:)); % This is the subset of data we will use in matrix form
[m,n] = size(rawdata); % We need m to know how many rows there are.
tempdata = [];

for ii = 1:(m-3); % Go down all the rows minus the ones that had text in rawdata.
    test = Finished(ii,1); % take the row
    if test == 1;  
        saveme = rawdata(ii+3,:); % Save the row if the test is passed. (If 0, fail.)
        tempdata = [tempdata; saveme]; % Concatenate new columns.
    end
end

% Engagement test.

% If participant is "Not engaged", then they get excluded.

% First step is to find all the engagement questions:

EngageTest = [];
for ii = 1:length(headers); % Go through every column in headers.
    pattern = "ngage"; % Look for this text pattern. There are three engagement questions and this finds them.
    TF = contains(headers{ii},pattern); % Logical... if you find the column, mark it.
    if TF == 1; % If true
        saveme = tempdata(:,ii); % Save the engagement column
       EngageTest = [saveme, EngageTest]; % And concatenate the engagement columns together.
    end    
end
 
EngageTest = cell2mat(EngageTest); % Converted to matrix to be easier.

[m,n] = size(tempdata); % m finds the row.
progressdata = [];

exclude = ismember(EngageTest,1); % If a person is 1 (not engaged), flag them as 1. All other values pass as 0.
for ii = 1:m; % Go through all the rows in tempdata
    test = exclude(ii,:); % Look at a row
    if test == 0; % Test to see if the row has any flags. If no flags, then good.
        saveme = tempdata(ii,:); % Save unflagged rows.
        progressdata = [progressdata; saveme];
    end    
end

% Comprehension test

% If you do not understand DG/UG Optimal Question, you get excluded.

% If participant answers greater than .5, they get excluded.

% First step is to find all Optimal Questions

OptimalTest = [];
for ii = 1:length(headers); % Go through every column in headers.
    pattern = "Optimal"; % Look for this text pattern. There are three engagement questions and this finds them.
    TF = contains(headers{ii},pattern); % Logical... if you find the column, mark it.
    if TF == 1; % If true
       saveme = progressdata(:,ii); % Save the engagement column
       OptimalTest = [saveme, OptimalTest]; % And concatenate the engagement columns together.
    end    
end
 
OptimalTest = cell2mat(OptimalTest); % Converted to matrix to be easier.

[m,n] = size(progressdata); % m finds the row.
progressdata2 = [];

OptimalTest(OptimalTest > .5) = 1; % Flag all values above .5 as 1
for ii = 1:m; % Go through all the rows in tempdata
     test = OptimalTest(ii,:); % Look at a row
     if test < .6; % Test to see if the row has any flags. If no flags, then good.
         saveme = progressdata(ii,:); % Save unflagged rows.
         progressdata2 = [progressdata2; saveme];
     end    
end

% 5) If you answer the engagement questions at the end of the survey
% incorrectly, you get excluded.

EngageTest2 = [];
for ii = 1:length(headers); % Go through every column in headers.
    pattern = "XYZ"; % Look for this text pattern. There are three engagement questions and this finds them.
    TF = contains(headers{ii},pattern); % Logical... if you find the column, mark it.
    if TF == 1; % If true
        saveme = progressdata2(:,ii); % Save the engagement column
       EngageTest2 = [saveme, EngageTest2]; % And concatenate the engagement columns together.
    end    
end

Test1 = cell2mat(EngageTest2(:,1)); % Converted to matrix to be easier.
Test2 = EngageTest2(:,2); % This one stays in cell format and will be tested differently.
Test3 = cell2mat(EngageTest2(:,3));

[m,n] = size(progressdata2); % m finds the row.
finaldata = [];

 for ii = 1:m; % Go through all the rows in progressdata2
      test = Test1(ii,:); % Look at a row
      if test == 4; % Test to see if the row has any flags. If no flags, then good.
          saveme = progressdata2(ii,:); % Save unflagged rows.
          finaldata = [finaldata; saveme];
      end    
 end
 
 finaldata2 = [];
 [m,n] = size(finaldata);
 
 for ii = 1:m; % Go through all the rows in finaldata
      test = Test3(ii,:); % Look at a row
      if test ~2; % Test to see if the row has any flags. If no flags, then good.
          saveme = finaldata(ii,:); % Save unflagged rows.
          finaldata2 = [finaldata2; saveme];
      end
 end

[m,n] = size(finaldata2);

indices = find(strcmp('1,2,4',Test2(:,1))); % This indexes the right answer from Test2
Subject = []; 
FinalData3 = [];
for ii = 1:length(indices)
    Subject = indices(ii,1); % Subject row.
    saveme = finaldata2(Subject,:); % Save the subject row.
    FinalData3 = [FinalData3; saveme]; % Concatenate the subjects 
end    

%%

data = [headers; FinalData3]; % Add the headers back so we know what we are working with.

% We now have the data that has been sorted through the exclusion criteria.

[m,n] = size(data);
TotalSubjects = m-1; % This sets the total subjects for further analysis.

% I left TotalSubjects in for each scale just in case subjects were
% arbitarily excluded for certain items. 

%% Machiavellian Scale

% Exclude all data except Machiavellian scores.

start = find(strcmp('Mach_questions_1',data(1,:))); % Find the Mach 1 column,
finish = find(strcmp('Mach_questions_20',data(1,:))); % Find the Mach 20 column.
N = 20; % Number of questions
IndexedColumns = round(linspace(start,finish, N)); % Index all of the Mach columns.
Mach_data = data(:,IndexedColumns); % Save them
Mach_data = cell2mat(Mach_data(2:end,:)); % This is the subset of Mach data we will use.

% Machiavellian Tactics Subscale

% Goal is to find a subject's tactical score

TacticScore = [];
TotalSubjects = TotalSubjects;
max_mach = 7; % This is for the reverse coding.
min_mach = 1; % This is for the reverse coding.
add_mach = max_mach + min_mach; 

for ii = 1:TotalSubjects % ii is the subject
    a = Mach_data(ii,1); 
    b = add_mach+(-1*(Mach_data(ii,4)));% Reverse code
    c = Mach_data(ii,6);
    d = Mach_data(ii,10);
    e = add_mach+(-1*(Mach_data(ii,12)));% Reverse code
    f = add_mach+(-1*(Mach_data(ii,15)));% Reverse code
    g = add_mach+(-1*(Mach_data(ii,16)));% Reverse code
    h = Mach_data(ii,18);
    i = add_mach+(-1*(Mach_data(ii,19))); % Reverse code
    total = a+b+c+d+e+f+g+h+i;
    TacticScore = [TacticScore, total];
end   

% Machiavellian Morality Subscale

MoralityScore= [];

for ii = 1:TotalSubjects % ii is the subject
    a = Mach_data(ii,8); 
    b = add_mach+(-1*(Mach_data(ii,20))); % Reverse code
    total = a+b;
    MoralityScore = [MoralityScore, total];
end   

% Machiavellian Views Subscale

ViewsScore= [];

for ii = 1:TotalSubjects; % ii is the subject
    a = add_mach+(-1*(Mach_data(ii,2)));% Reverse code
    b = Mach_data(ii,3);
    c = Mach_data(ii,5);
    d = Mach_data(ii,7);
    e = add_mach+(-1*(Mach_data(ii,9)));% Reverse code
    f = Mach_data(ii,11);
    g = add_mach+(-1*(Mach_data(ii,13)));% Reverse code
    h = add_mach+(-1*(Mach_data(ii,14)));% Reverse code
    i = Mach_data(ii,17);
    total = b+a+c+d+e+f+g+h+i;
    ViewsScore = [ViewsScore, total];
end   

Overall_Mach_Score = ((TacticScore + MoralityScore + ViewsScore)); %Overall Mach score is the mean of the 20 items (after reversing the appropriate items).

%% Personal Norm Scale

% Find the columns you will need.

start = find(strcmp('PNR_Scale_1',data(1,:))); % Process is same as in Mach.
finish = find(strcmp('PNR_Scale_9',data(1,:)));
N = 9; % Number of questions
IndexedColumns = round(linspace(start,finish, N));
PNR_data = data(:,IndexedColumns);
PNR_data = cell2mat(PNR_data(2:end,:)); % This is the subset of data we will use.

% PNR Scale

PNRScore = [];
TotalSubjects = TotalSubjects;
max_pnr = 7;
min_pnr = 1;
add_pnr = max_pnr + min_pnr;

for ii = 1:TotalSubjects % ii is the subject
    a = PNR_data(ii,1); 
    b = PNR_data(ii,2);
    c = PNR_data(ii,3);
    d = PNR_data(ii,4);
    e = PNR_data(ii,5);
    f = PNR_data(ii,6);
    g = PNR_data(ii,7);
    h = PNR_data(ii,8);
    i = PNR_data(ii,9);
    total = a+b+c+d+e+f+g+h+i;
    PNRScore = [PNRScore, total];
end   

%% EQ Scale

% Find the columns you will need.

start = find(strcmp('EQ_1',data(1,:)));
finish = find(strcmp('EQ_30',data(1,:)));
N = 30; % Number of questions
IndexedColumns = round(linspace(start,finish, N));
EQ_data = data(:,IndexedColumns);
EQ_data = cell2mat(EQ_data(2:end,:)); % This is the subset of data we will use.

WellbeingScore = [];
SelfcontrolScore = [];
EmotionalityScore = [];
SociabilityScore = [];
TotalEQScore = [];
TotalSubjects = TotalSubjects; % Put in total subjects
max_EQ = 7;
min_EQ = 1;
add_EQ = max_EQ + min_EQ;

for ii = 1:TotalSubjects % ii is the subject
    EQ1 = EQ_data(ii,1);
    EQ2 = add_EQ+(-1*(EQ_data(ii,2)));% Reverse code
    EQ3 = EQ_data(ii,3);
    EQ4 = add_EQ+(-1*(EQ_data(ii,4)));% Reverse code
    EQ5 = add_EQ+(-1*(EQ_data(ii,5)));% Reverse code
    EQ6 = EQ_data(ii,6);
    EQ7 = add_EQ+(-1*(EQ_data(ii,7)));% Reverse code
    EQ8 = add_EQ+(-1*(EQ_data(ii,8)));% Reverse code
    EQ9 = EQ_data(ii,9);
    EQ10 = add_EQ+(-1*(EQ_data(ii,10)));% Reverse code
    EQ11 = EQ_data(ii,11);
    EQ12 = add_EQ+(-1*(EQ_data(ii,12)));% Reverse code
    EQ13 = add_EQ+(-1*(EQ_data(ii,13)));% Reverse code
    EQ14 = add_EQ+(-1*(EQ_data(ii,14)));% Reverse code
    EQ15 = EQ_data(ii,15);
    EQ16 = add_EQ+(-1*(EQ_data(ii,16)));% Reverse code
    EQ17 = EQ_data(ii,17);
    EQ18 = add_EQ+(-1*(EQ_data(ii,18)));% Reverse code
    EQ19 = EQ_data(ii,19);
    EQ20 = EQ_data(ii,20);
    EQ21 = EQ_data(ii,21);
    EQ22 = add_EQ+(-1*(EQ_data(ii,22)));% Reverse code 
    EQ23 = EQ_data(ii,23); 
    EQ24 = EQ_data(ii,24);
    EQ25 = add_EQ+(-1*(EQ_data(ii,25)));% Reverse code
    EQ26 = add_EQ+(-1*(EQ_data(ii,26)));% Reverse code
    EQ27 = EQ_data(ii,27);
    EQ28 = add_EQ+(-1*(EQ_data(ii,28)));% Reverse code
    EQ29 = EQ_data(ii,29); 
    EQ30 = EQ_data(ii,30); 
    
    total = EQ1+EQ2+EQ3+EQ4+EQ5+EQ6+EQ7+EQ8+EQ9+EQ10+EQ11+EQ12+EQ13+EQ14+EQ15+EQ16+EQ17+EQ18+EQ19+EQ20+EQ21+EQ22+EQ23+EQ24+EQ25+EQ26+EQ27+EQ28+EQ29+EQ30;
    TotalEQScore = [TotalEQScore, total];
    
    % Subscales are calculated from Deshawn's code

    WS = EQ5 + EQ20 + EQ9 + EQ24 + EQ12 + EQ27;
    WellbeingScore = [WellbeingScore, WS]; 
    
    SC = EQ4 + EQ19 + EQ7 +EQ22 + EQ15 + EQ30;
    SelfcontrolScore = [SelfcontrolScore, SC];
    
    ES = EQ1 + EQ16 + EQ2 + EQ17 + EQ8 + EQ23 + EQ13 + EQ28;
    EmotionalityScore = [EmotionalityScore, ES];
    
    SS = EQ6 + EQ21 + EQ10 + EQ25 + EQ11 + EQ26;
    SociabilityScore = [SociabilityScore, SS]; 
end 

%% Demographics

% Find the columns you will need.

start = find(strcmp('Age',data(1,:))); % Take slice of data.
finish = find(strcmp('Education',data(1,:)));
N = 6; % Number of questions
IndexedColumns = round(linspace(start,finish, N));
demo_data = data(:,IndexedColumns);
demo_data = cell2mat(demo_data(2:end,:)); % This is the subset of data we will use.

age = demo_data(:,1);

gender = demo_data(:,2);

males = sum(gender(:) == 1);
females = TotalSubjects - males;

race = demo_data(:,4);

white = sum(race(:) == 1);
black = sum(race(:) == 2);
hispanic = sum(race(:) == 3);
indian = sum(race(:) == 4);
asian = sum(race(:) == 5);
PI = sum(race(:) == 6);
mixed = sum(race(:) == 7);
other = (race - white - black - hispanic - indian - asian - PI - mixed);

education = demo_data(:,6);

middle_school = sum(education(:) == 1);
high_school = sum(education(:) == 2);
training = sum(education(:) == 3);
associate = sum(education(:) == 4);
bachelor = sum(education(:) == 5);
master = sum(education(:) == 6);
doctor = sum(education(:) == 7); % No way
professional = sum(education(:) == 8);

%% Clean UG-A Data 

UG_A_temp = [];

for ii = 1:length(headers); % Go through all the names in the header.
    pattern = "UG_A_"; % Find strings with this information.
    TF = contains(headers{ii},pattern); % True or false, does the string have UG_A in it?
    if TF == 1; % if it does...
        saveme= data(:,ii); % Save the column...
        UG_A_temp = [saveme, UG_A_temp]; % In UG_A_Temp.
    end    
end     
    
UG_A = cell2mat(UG_A_temp(2:end,:)); % Change temp into a matrix for analysis.
    
UG_A_Mean_Offers = mean(UG_A'); % Find the mean offers for all the participants.

%% Clean DG-A Data 

DG_A_temp = [];

for ii = 1:length(headers); % See process in UG-A for comments.
    pattern = "DG_A_";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        DG_A_temp = [saveme, DG_A_temp];
    end    
end     
    
DG_A = cell2mat(DG_A_temp(2:end,:));
    
DG_A_Mean_Offers = mean(DG_A');

%% Clean UG-B Data 

UG_B_temp = [];

for ii = 1:length(headers); % See process in UG-A for comments. 
    pattern = "UG_B_X";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        UG_B_temp = [saveme, UG_B_temp];
    end    
end     
    
UG_B = cell2mat(UG_B_temp(2:end,:)); % Just in case, I have the raw data in all the columns here.

% 2 is accept, 1 is reject % MAKE SURE THIS IS RIGHT!!!!!!

% Let's group all the offers under their respective offer amount.

X6 = [];
X5 = [];
X4 = []; 
X3 = [];
X2 = [];
X1 = [];
X0 = [];

for ii = 1:length(headers); % See process in UG-A
    pattern = "UG_B_X6";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        X6 = [saveme, X6];
    end    
end   

for ii = 1:length(headers); % ditto
    pattern = "UG_B_X5";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        X5 = [saveme, X5];
    end    
end 

for ii = 1:length(headers); % ditto
    pattern = "UG_B_X4";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        X4 = [saveme, X4];
    end    
end 
for ii = 1:length(headers); % ditto
    pattern = "UG_B_X3";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        X3 = [saveme, X3];
    end    
end 

for ii = 1:length(headers); % ditto
    pattern = "UG_B_X2";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        X2 = [saveme, X2];
    end    
end 

for ii = 1:length(headers); % ditto
    pattern = "UG_B_X1";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        X1 = [saveme, X1];
    end    
end 

for ii = 1:length(headers); % ditto
    pattern = "UG_B_X0";
    TF = contains(headers{ii},pattern);
    if TF == 1;
        saveme = data(:,ii);
        X0 = [saveme, X0];
    end    
end 

%% Behavioral Survey answers

UG_Fair_Column = find(strcmp('UG_Fair_1',data(1,:))); % Find the UG_Fair string in the data file.
UG_Fair = data(:,UG_Fair_Column); % Save it.

UG_Optimal_Column = find(strcmp('UG_Optimal_1',data(1,:))); % etc.
UG_Optimal = data(:,UG_Optimal_Column);

UG_Ask_Column = find(strcmp('UG_Ask',data(1,:)));
UG_Ask = data(:,UG_Ask_Column);

DG_Fair_Column = find(strcmp('DG_Fair_1',data(1,:)));
DG_Fair = data(:,DG_Fair_Column);

DG_Optimal_Column = find(strcmp('DG_Optimal_1',data(1,:)));
DG_Optimal = data(:,DG_Optimal_Column);

DG_Ask_Column = find(strcmp('DG_Ask',data(1,:)));
DG_Ask = data(:,DG_Ask_Column);

%% Save Desired Variables

save('Behavioral_Data', 'UG_A', 'DG_A', 'UG_A_Mean_Offers','DG_A_Mean_Offers', 'X0', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6')
save('Behavioral_Survey_Data', 'DG_Optimal', 'DG_Fair', 'DG_Ask', 'UG_Optimal', 'UG_Fair', 'UG_Ask')
save('Survey_Data', 'PNRScore', 'Overall_Mach_Score', 'ViewsScore', 'MoralityScore', 'TacticScore', 'TotalEQScore', 'WellbeingScore', 'SelfcontrolScore', 'EmotionalityScore', 'SociabilityScore')
save('Demographic_Data', 'age', 'gender', 'males', 'females', 'race', 'white', 'black', 'hispanic', 'indian', 'asian', 'PI', 'mixed', 'other', 'education', 'middle_school', 'high_school', 'training', 'associate', 'bachelor', 'master', 'doctor', 'professional')