%% Initialization

clear all;
close all;
clc;

[n,t,rawdata] = xlsread('Exponentials_1117.csv'); % This converts the CSV file into a matlab file.

% Rawdata needs to be processed through the exclusion criteria. Then it
% becomes "data"

%% Briefing

% This script will clean up the data from the Stay/Go Qualtrics survey

% Daniel Sazhin
% 07-25-2020
% DVS Lab
% Temple University

% This script takes the CSV file and cleans the data. At the end of it, you
% will have the following useful information from our participants in
% several files:

% INSERT VARIABLES HERE

%% Headers

use = t(1,:); % Pulls out the headers
headers = cell2table(use); % This turns the headers into a table format
headers = table2array(headers);

%% Exclusion criteria

% Commitment to quality questionnaire. If rejected, excluded.
% Fail to finish, excluded ???? not pre-registered
% Random attention checks. If failed, excluded.


%%%%%% Not commpleted yet %%%%%%%%%%%
% Exclude subjects who complete the survey in (>3SD) of the median completion time.
% Exclude subjects who leave at 0 turns > 50 trials (do not understand the task)
% Exclude subjects who leave at 10 turns > 50 trials (do not understand the task)


%% Fail to finish, exclude

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


%% Attention Check 1

tempdata2 = [];
temp = find(strcmp('Attention Check 1',rawdata(1,:))); % This finds the Finished column.
Finished = tempdata(:,temp); % Takes out the column from the data table.
Finished = cell2mat(Finished(4:end,:)); % This is the subset of data we will use in matrix form
[m,n] = size(tempdata); % We need m to know how many rows there are.


for ii = 1:(m-3); % Go down all the rows minus the ones that had text in rawdata.
    test = Finished(ii,1); % take the row
    if test == 7;
        saveme = tempdata(ii+3,:); % Save the row if the test is passed. (If 0, fail.)
        tempdata2 = [tempdata2; saveme]; % Concatenate new columns.
    end
end

%% Attention Check 2

tempdata3 = [];
temp = find(strcmp('Attention Check 2',rawdata(1,:))); % This finds the Finished column.
Finished = tempdata2(:,temp); % Takes out the column from the data table.
Finished = cell2mat(Finished(4:end,:)); % This is the subset of data we will use in matrix form
[m,n] = size(tempdata2); % We need m to know how many rows there are.


for ii = 1:(m-3); % Go down all the rows minus the ones that had text in rawdata.
    test = Finished(ii,1); % take the row
    if test == 8;
        saveme = tempdata2(ii+3,:); % Save the row if the test is passed. (If 0, fail.)
        tempdata3 = [tempdata3; saveme]; % Concatenate new columns.
    end
end

%% Attention Check 2

tempdata4 = [];
temp = find(strcmp('Duration (in seconds)',rawdata(1,:))); % This finds the Duration column.
Duration = tempdata3(:,temp); % Takes out the column from the data table.
Duration = cell2mat(Duration); % This is the subset of data we will use in matrix form
[m,n] = size(tempdata3); % We need m to know how many rows there are.
median_Duration = median(Duration);
std_Duration = std(Duration);
test = median_Duration + 3*std_Duration;
% Exclude subjects who complete the survey in (>3SD) of the median completion time

for ii = 1:(m); % Go down all the rows minus the ones that had text in rawdata.
    test_row = Duration(ii,1); % take the row
    if test_row < test;
        saveme = tempdata2(ii+3,:); % Save the row if the test is passed. (If 0, fail.)
        tempdata4 = [tempdata4; saveme]; % Concatenate new columns.
    end
end

%%

data = [headers; tempdata4]; % Add the headers back so we know what we are working with.

% We now have the data that has been sorted through the exclusion criteria.

[m,n] = size(data);
TotalSubjects = m-1; % This sets the total subjects for further analysis.

% I left TotalSubjects in for each scale just in case subjects were
% arbitarily excluded for certain items.

%% Risk Scale

% Exclude all data except risk aversion scores.

start = find(strcmp('Risk aversion_1',data(1,:))); % Find the risk 1 column,
finish = find(strcmp('Risk aversion_12',data(1,:))); % Find the risk 12 column.
N = 12; % Number of questions
IndexedColumns = round(linspace(start,finish, N)); % Index all of the Risk columns.
Risk_data = data(:,IndexedColumns); % Save them
Risk_data = cell2mat(Risk_data(2:end,:)); % This is the subset of Risk data we will use.

Normalized_Risk = Risk_data - 1;
save_risk = [];

for ii = 1:length(Normalized_Risk)
    sub = Normalized_Risk(ii,:);
    
    x = rand(12, 1);
    y = Normalized_Risk(ii,:);
    y(1:end) = x(1:end) > 0.3; % To avoid perfect seperation
    % Fit model
    mdl = fitglm(x, y, "Distribution", "binomial");
    xnew = linspace(0,1,1200)'; % test data
    ynew = predict(mdl, xnew);
    %scatter(x, y);
    %hold on;
    %plot(xnew, ynew);
    
    risk_logistic = min(find(ynew>.5));
    
    save_risk = [save_risk; risk_logistic/1200];

end

%% Loss aversion

% Exclude all data except loss aversion scores.

start = find(strcmp('Loss aversion_1',data(1,:))); % Find the Loss 1 column,
finish = find(strcmp('Loss aversion_12',data(1,:))); % Find the Loss 12 column.
N = 12; % Number of questions
IndexedColumns = round(linspace(start,finish, N)); % Index all of the Loss columns.
Loss_data = data(:,IndexedColumns); % Save them
Loss_data = cell2mat(Loss_data(2:end,:)); % This is the subset of Loss data we will use.

% Make a logistic regression to assess indifference

Normalized_Loss = Loss_data - 1;

save_loss = [];

for ii = 1:length(Normalized_Loss)
    sub = Normalized_Loss(ii,:)
    
    x = rand(12, 1);
    y = Normalized_Loss(ii,:);
    y(1:end) = x(1:end) > 0.3; % To avoid perfect seperation
    % Fit model
    mdl = fitglm(x, y, "Distribution", "binomial");
    xnew = linspace(0,1,1200)'; % test data
    ynew = predict(mdl, xnew);
    %scatter(x, y);
    %hold on;
    %plot(xnew, ynew);
    
    loss_logistic = min(find(ynew>.5));
    
    save_loss = [save_loss; loss_logistic/1200];

end





%% 7up7down

% Exclude all data except 7up7down scores.

start = find(strcmp('7up7down_1',data(1,:))); % Find the 7up7down 1 column,
finish = find(strcmp('7up7down_14',data(1,:))); % Find the 7up7down 12 column.
N = 14; % Number of questions
IndexedColumns = round(linspace(start,finish, N)); % Index all of the SevenUp columns.
SevenUpSevenDown_data = data(:,IndexedColumns); % Save them

SevenUpSevenDown_data = strrep(SevenUpSevenDown_data,'Never or hardly ever','1');
SevenUpSevenDown_data = strrep(SevenUpSevenDown_data,'Sometimes','2');
SevenUpSevenDown_data = strrep(SevenUpSevenDown_data,'Often','3');
SevenUpSevenDown_data = strrep(SevenUpSevenDown_data,'Very often or almost constantly','4');


[n,m] = size(SevenUpSevenDown_data);
SevenUpSevenDown_data_process = [];

for ii = 1:(n-1) % for each participant
    participant = [];
    for jj = 1:m % for each question
        response = SevenUpSevenDown_data(ii+1,jj);
        response = str2double(response);
        participant = [participant, response];
    end
    SevenUpSevenDown_data_process = [SevenUpSevenDown_data_process; participant]; % Final scores
end

SevenUp = [];
[n,m] = size(SevenUpSevenDown_data_process);
for ii = 1:n
    row = SevenUpSevenDown_data_process(ii,:);
    participant = row([1 3 4 6 7 8 13]);
    saveme = sum(participant);
    SevenUp = [SevenUp; saveme];
end

SevenDown = [];
[n,m] = size(SevenUpSevenDown_data_process);
for ii = 1:n
    row = SevenUpSevenDown_data_process(ii,:);
    participant = row([2 5 9 10 11 12 14]);
    saveme = sum(participant);
    SevenDown = [SevenDown; saveme];
end

%% AADIS

% Note: This data might be less stable. There was quite a bit of
% hard-coding with respect to recoding text into integers... other
% questions might not have been coded in correctly.

% Exclude all data except AADIS

start = find(strcmp('Q1',data(1,:))); % Find the Loss 1 column,
finish = find(strcmp('Q14',data(1,:))); % Find the Loss 12 column.
N = 14; % Number of questions
IndexedColumns = round(linspace(start,finish, N)); % Index all of the AADIS columns.
AADIS_data = data(:,IndexedColumns); % Save them

% Change Nan to 0 (Non-answered questions)

[n,m] = size(AADIS_data);
for ii = 1:n % for every participant
    for jj = 1:m % find a cell
    test = cell2mat(AADIS_data(ii,jj));
    if isnan(test) == 1
        AADIS_data(ii,jj) = {'0'};
    end
    end
end
        

% Question 1

AADIS_data = strrep(AADIS_data,'never','0');
AADIS_data = strrep(AADIS_data,'once or twice a year','2');
AADIS_data = strrep(AADIS_data,'once or twice a  month','3');
AADIS_data = strrep(AADIS_data,'every weekend','4');
AADIS_data = strrep(AADIS_data,'several times a week','5');
AADIS_data = strrep(AADIS_data,'every day','6');
AADIS_data = strrep(AADIS_data,'several times a day','7');

% Question 2

AADIS_data = strrep(AADIS_data,'used alcohol or drugs','0');
AADIS_data = strrep(AADIS_data,'not for over a year','2');
AADIS_data = strrep(AADIS_data,'between 6  months and 1 year [before]','3');
AADIS_data = strrep(AADIS_data,'several weeks ago [before]','4');
AADIS_data = strrep(AADIS_data,'last week [the week before]','5');
AADIS_data = strrep(AADIS_data,'yesterday [the day before]','6');
AADIS_data = strrep(AADIS_data,'today','7');

% Question 3

AADIS_data = strrep(AADIS_data,'I like the feeling','1');
AADIS_data = strrep(AADIS_data,'to be like my friends','2');
AADIS_data = strrep(AADIS_data,'I am bored; or just to have fun','3');
AADIS_data = strrep(AADIS_data,"kickin' it",'3');
AADIS_data = strrep(AADIS_data,'3 ("3")','3');
AADIS_data = strrep(AADIS_data,'I feel stressed, nervous, tense, full of worries or problems','4');
AADIS_data = strrep(AADIS_data,'4,I feel sad, lonely, sorry for myself','4');
AADIS_data = strrep(AADIS_data,'I am bored; or just to have fun','5');

% Question 4

AADIS_data = strrep(AADIS_data,'wine','1');
AADIS_data = strrep(AADIS_data,'beer','2');
AADIS_data = strrep(AADIS_data,'mixed drinks','3');
AADIS_data = strrep(AADIS_data,'hard liquor (vodka, whisky, etc.)','4');
AADIS_data = strrep(AADIS_data,'a substitute for alcohol','5');
AADIS_data = strrep(AADIS_data,'a  substitute for alcohol','5');

% Question 5

AADIS_data = strrep(AADIS_data,'Supervised by parents or relatives','1');
AADIS_data = strrep(AADIS_data,'from brothers or sisters','2');
AADIS_data = strrep(AADIS_data,"from home without parents' knowledge",'3');
AADIS_data = strrep(AADIS_data,'get from friends','4');
AADIS_data = strrep(AADIS_data,'buy my own (on the street or with false 10)','5');

% Question 6

AADIS_data = strrep(AADIS_data,'never','0');
AADIS_data = strrep(AADIS_data,'after age 15','2');
AADIS_data = strrep(AADIS_data,'at ages 14 or 15','3');
AADIS_data = strrep(AADIS_data,'at ages 12 or 13','4');
AADIS_data = strrep(AADIS_data,'at ages 10 or 11','5');
AADIS_data = strrep(AADIS_data,'before age 10','6');

% Question 7

AADIS_data = strrep(AADIS_data,'at night','1');
AADIS_data = strrep(AADIS_data,'afternoons/after school','2');
AADIS_data = strrep(AADIS_data,'before or during school or work','3');
AADIS_data = strrep(AADIS_data,'in the morning or when I first awaken','4');
AADIS_data = strrep(AADIS_data,'I often get up during my sleep to use alcohol or drugs','5');

% Question 8

AADIS_data = strrep(AADIS_data,'curiosity','1');
AADIS_data = strrep(AADIS_data,'parents or relatives offered','2');
AADIS_data = strrep(AADIS_data,'friends encouraged me; to have fun','3');
AADIS_data = strrep(AADIS_data,'to get away from my problems','4');
AADIS_data = strrep(AADIS_data,'to get high or drunk','5');

% Question 9

AADIS_data = strrep(AADIS_data,'1 drink','1');
AADIS_data = strrep(AADIS_data,'2 drinks','2');
AADIS_data = strrep(AADIS_data,'3-4 drinks','3');
AADIS_data = strrep(AADIS_data,'5-9 drinks','4');
AADIS_data = strrep(AADIS_data,'10 or more drinks','5');

% Question 10

AADIS_data = strrep(AADIS_data,'parents or adult relatives','1');
AADIS_data = strrep(AADIS_data,'with brothers or sisters','2');
AADIS_data = strrep(AADIS_data,'with friends or relatives own age','3');
AADIS_data = strrep(AADIS_data,'with older friends','4');
AADIS_data = strrep(AADIS_data,'alone','5');

% Question 11

AADIS_data = strrep(AADIS_data,'loose, easy feeling','1');
AADIS_data = strrep(AADIS_data,'got moderately high','2');
AADIS_data = strrep(AADIS_data,'got drunk or wasted','3');
AADIS_data = strrep(AADIS_data,'became ill','4');
AADIS_data = strrep(AADIS_data,'passed out or overdosed','5');
AADIS_data = strrep(AADIS_data,"used a lot and next day didn't remember what happened",'6');

% Question 12

AADIS_data = strrep(AADIS_data,'none','0');
AADIS_data = strrep(AADIS_data,'has interfered with talking to someone','2');
AADIS_data = strrep(AADIS_data,'has prevented me from having a good time','3');
AADIS_data = strrep(AADIS_data,'has interfered with my school work','4');
AADIS_data = strrep(AADIS_data,'have lost friends because of use','5');
AADIS_data = strrep(AADIS_data,'has gotten me into trouble at home','6');
AADIS_data = strrep(AADIS_data,'was in a fight or destroyed property','7');
AADIS_data = strrep(AADIS_data,'has resulted in an accident, an arrest, or being punished at school for using alcohol or drugs','8');

% Question 13

AADIS_data = strrep(AADIS_data,'no problem at all','0');
AADIS_data = strrep(AADIS_data,'I can control it and set limits on myself','1');
AADIS_data = strrep(AADIS_data,'I can control myself, but my friends easily influence me','3');
AADIS_data = strrep(AADIS_data,'I often feel bad about my use','4');
AADIS_data = strrep(AADIS_data,'I need help to control myself','5');
AADIS_data = strrep(AADIS_data,'I have had professional help to control my drinking or drug use','6');

% Question 14
AADIS_data = strrep(AADIS_data,'say or normal for my age','0');
AADIS_data = strrep(AADIS_data,'canâ€™t 0','0');
AADIS_data = strrep(AADIS_data,'when I use I tend to neglect my family or friends','2');
AADIS_data = strrep(AADIS_data,'my family or friends advise me to control or cut down on my use','3');
AADIS_data = strrep(AADIS_data,'my family or friends tell me to get help for my alcohol or drug use','4');
AADIS_data = strrep(AADIS_data,'my family or friends have already gone for help about my alcohol or drug use','5');

[n,m] = size(AADIS_data);
AADIS_data_final = [];

for ii = 1:(n-1) % for each participant
    participant = [];
    for jj = 1:m % for each question
        response = AADIS_data(ii+1,jj);
        
        response = str2num(response);
        participant = [participant, response];
        participant_sum = sum(participant);
    end
    AADIS_data_final = [AADIS_data_final; participant_sum]; % Final scores
end


%% Intolerance of Uncertainty

% Exclude all data except loss aversion scores.

start = find(strcmp('Uncertainty_1',data(1,:))); % Find the Loss 1 column,
finish = find(strcmp('Uncertainty_12',data(1,:))); % Find the Loss 12 column.
N = 12; % Number of questions
IndexedColumns = round(linspace(start,finish, N)); % Index all of the Loss columns.
Uncertainty_data = data(:,IndexedColumns); % Save them

% 1) Find and replace text with values.

Uncertainty_data = strrep(Uncertainty_data,'Not at all characteristic of me','1');
Uncertainty_data = strrep(Uncertainty_data,'A little characteristic of me','2');
Uncertainty_data = strrep(Uncertainty_data,'Somewhat characteristic of me','3');
Uncertainty_data = strrep(Uncertainty_data,'Very characteristic of me','4');
Uncertainty_data = strrep(Uncertainty_data,'Entirely characteristic of me','5');

[n,m] = size(Uncertainty_data);
Uncertainty_data_final = [];

for ii = 1:(n-1) % for each participant
    participant = [];
    for jj = 1:m % for each question
        response = Uncertainty_data(ii+1,jj);
        response = cell2mat(response);
        response = str2num(response);
        participant = [participant, response];
    end
    Uncertainty_data_final = [Uncertainty_data_final; participant];
end

Uncertainty_data_final = sum(Uncertainty_data_final,2);


%% Demographics

% Find the columns you will need.

start = find(strcmp('First',data(1,:))); % Take slice of data.
finish = find(strcmp('Income',data(1,:)));
N = 9; % Number of questions
IndexedColumns = round(linspace(start,finish, N));
demo_data = data(:,IndexedColumns);
Age = demo_data(2:end,3);
Age=cell2mat(Age);
mean_age= mean(Age);
std_age= std(Age);
range_age=range(Age);

Ethnicity = demo_data(2:end,7);
Gender = demo_data(2:end,5);
[N,~] = size(find(strcmp('Female',Gender)));
[M,~] = size(find(strcmp('Hispanic or Latino',Ethnicity)));
[O,~] = size(Ethnicity);
prop_hispanic = M/O;

%% Behavioral data

% This section is a little bit complicated. There are two ways that a
% participant could have went through a trial... they could have Left or
% Continued through the 10th trial.

% To figure out their decisions, we need to select the trials they Left and
% those respective indices. And the trials they Continued.

% First step, find Q9 indices. This will let us test whether a subject
% stayed through 10th turn later.

indices_Q9 =[];
k = strfind(data(1,:),'_Q9_','ForceCellOutput',true);
for ii = 1:length(k)
    value = cell2mat(k(ii));
    saverow = [];
    if value > 0
        saverow = ii;
    end
    indices_Q9 = [indices_Q9,saverow];
end



% Now we will find the trial behaviors:


for ii = 1:length(data(2:end,1)) % For each participant
    
    indices = find(strcmp('Leave?',data(ii+1,:)))'; % Take their indices when they Left
    Participant_Matrix = [];
    Participant = [];
    
    for jj = 1:length(indices) % Pick a trial
        index = indices(jj);
        saveme = cell2mat(data(1,index)); % This saves the header.
        
        % Decode saveme %
        
        Alpha = str2double(extractBetween(saveme,"D","_E"));
        Final_Value = str2double(extractBetween(saveme,"_E","_Q"));
        Turn_Left = str2double(extractBetween(saveme,"_Q","_"));
        Trial_Earning = extractAfter(saveme,"_Q");
        Trial_Earning = str2double(extractAfter(Trial_Earning,"_"));
        
        if Turn_Left == 10 % Determines trial earnings
            Trial_Earning = Final_Value;
        else
            Trial_Earning = 10 - Trial_Earning;
        end
        
        Trial = [Alpha,Final_Value,Turn_Left,Trial_Earning]; % Saves a trial as raw info
        Participant_Matrix = [Participant_Matrix; Trial];
        Participant = array2table(Participant_Matrix(1:end,:),'VariableNames', {'Alpha' 'Final_Value' 'Turn_Left' 'Trial_earnings'});
    end

    % Now we add in turn 10 decisions (staying all the way through)
    
    for kk = 1:length(indices_Q9) % Find all Turn 9 indices.
        index = indices_Q9(kk); % Save that index for testing.
        saveme = cell2mat(data(1,index)); % Save the header.
        
        index_test= indices_Q9(kk); %
        test = find(strcmp('Continue to next turn?',data(ii+1,index_test)))'; % Did the participant Continue on Q9??
        
        if test > 0 % If so...
            saveme = cell2mat(data(1,index)); % Save the header
            
            % Decode the header
            
            Alpha = str2double(extractBetween(saveme,"D","_E"));
            Final_Value = str2double(extractBetween(saveme,"_E","_Q"));
            Turn_Left = 10;
            Trial_Earning = extractAfter(saveme,"_Q");
            Trial_Earning = str2double(extractAfter(Trial_Earning,"_"));
            
            Trial_Earning = Final_Value;
            
            Trial = [Alpha,Final_Value,Turn_Left,Trial_Earning]; % Combine nicely
            Participant_Matrix = [Participant_Matrix; Trial];
            Participant = array2table(Participant_Matrix(1:end,:),'VariableNames', {'Alpha' 'Final_Value' 'Turn_Left' 'Trial_earnings'});
            
            
        end
    end
    
    filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
    writetable(Participant, filename) % Save file in table form
    
    test_high = (sum(Participant.Turn_Left) > 5000); % Test if participant leaves 50 turns at turn 10
    test_low = (sum(Participant.Turn_Left) < 60); % Test if participant leaves 50 turns at turn 1
    
        if test_high == 1 % print participant number
           disp(ii)
           disp('failed test high')
        end
    
        if test_low == 1 % print participant number
            disp(ii)
            disp('failed test low')
        end
    
end


  
    
    
    



%% Save

save('Survey_Data', 'Uncertainty_data_final', 'AADIS_data_final', 'SevenDown', 'SevenUp', 'save_risk', 'save_loss');
