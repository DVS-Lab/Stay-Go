clear all
close all
clc

%% Cleanup 2

load data_final.mat

data = data_final;

%% Rerun behavioral data to order it nicely for analysis

% This section is a little bit complicated. There are two ways that a
% participant could have went through a trial... they could have Left at any time or
% Continued through the 10th trial.

% To figure out their decisions, we need to select the trials they Left and
% those respective indices. And the trials they Continued.

% First step, find Q9 indices. This will let us test whether a subject
% stayed through 10th turn.

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

Participant = [];
testthis = [];
save_exclusion = [];

for ii = 1:length(data(2:end,1))'; % For each participant
    
    indices = find(strcmp('Leave?',data(ii+1,:)))'; % Take their indices when they Left
    Participant_Matrix = [];
    savematrix = [];
    for jj = 1:length(indices); % Pick a trial
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
        
    end
    
    for yy = 1:length(Participant_Matrix)
        row = Participant_Matrix(yy,3);
        if row > 1
            savematrix = [savematrix; Participant_Matrix(yy,:)];
        end
    end
    
    Participant = array2table(savematrix(1:end,:),'VariableNames', {'Alpha' 'Final_Value' 'Turn_Left' 'Trial_earnings'});
 
    filename = ['Participant_Matrix_' sprintf('%01d',ii) '.csv'];
    writetable(Participant, filename) % Save file in table form
end


  
%% Save

save('Survey_Data_robust')
