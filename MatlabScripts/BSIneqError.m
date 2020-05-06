% Reads in data
format long;
value = csvread('BSIneqError.csv');

% Retrieves h and error values
h = value(1:end,1);
error = value(1:end,2);

% Create plot
figure(1);
loglog(h,error,'m','LineWidth',2);
xlabel('Step size h'); ylabel('Error norm');
title('Maximum error norm for American put option value against spatial step-size'); % Improve with t
ax = gca; ax.FontSize = 14; axis tight;


