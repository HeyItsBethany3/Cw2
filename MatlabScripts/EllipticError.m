% Reads in data
format long;
value = csvread('EllipticError.csv');

% Retrieves h and error values
h = value(1:end,1);
error = value(1:end,2);

% Create plot
figure(1);
loglog(h,error,'m','LineWidth',2);
xlabel('log(h)'); ylabel('log(error)');
title('Grid error norm for Elliptic PDE against spatial step-size h');
ax = gca; ax.FontSize = 14; axis tight;

% Displays gradient 
fitlm(log(h), log(error))
