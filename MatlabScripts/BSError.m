% Reads in data
format long;
value = csvread('BSError.csv');

% Retrieves h and error values
h = value(1:end,1);
error = value(1:end,2);

% Create plot
figure(1);
loglog(h,error,'m','LineWidth',2);
xlabel('log(h)'); ylabel('log(error)');
title('Maximum error norm for European put option value against step-size h'); % Improve with t
ax = gca; ax.FontSize = 14; axis tight;

% Displays gradient 
fitlm(log(h), log(error))

