clear all 
close all
clc

%% Distributions

% Selected parameters

X0= 0; % Start value
N = 10; % Number of turns 
Final_Values = [5:15]; % Cost at final turn
alphas = [1.4, 1.7, 2.0, 2.3, 2.6]; % Growth Factors

Distributions = []; % Preallocate

for aa = 1:length(alphas) % Choose growth factor

for jj = 1:length(Final_Values) % Choose final value

Xn= Final_Values(jj); % Final value
alpha = alphas(aa); % Growth Factor 

denom = [];
for ii = 1:(N-1)
    saveme = (alpha^ii);
    denom = [denom, saveme]; % The equation has a recursive element...
end

m = (Xn - (X0*(alpha^N)))/(sum(denom)+1); % m is the linear normalizing factor which constrains the function to the starting point 0.

Exponential_Curve = [];

for ii = 1:N;
Xk = (X0)*(alpha^ii) + m*(sum(denom(1:(ii-1)))+1); % This is the guiding equation. Xk selects a turn
Exponential_Curve = [Exponential_Curve, Xk];
end

Distributions = [Distributions; Exponential_Curve]; 

end
end


%% Figures


figure

for kk = 1:length(Distributions)
  
    X = [1:N];
    plot(X,Distributions(kk,:));
    xlabel('Turns')
    ylabel('Dollars')
    title('Exponential Distribution')
    ylim([0 Final_Values(jj)])
    
    hold on
end

%% Export

name = 'Exponentials.csv';
writematrix(Distributions, name); % Save as csv file
