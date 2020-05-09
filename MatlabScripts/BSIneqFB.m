% Reads in data
format long;
value = csvread('BSIneqFB.csv');

x = value(1,1:end-1);  % x values
fb = value(2, 1:end-1); % free boundary t values

% Create plots
figure(1);
plot(fb, x,'m');
title('Optimal stopping boundary for each node x');
xlabel('Time to maturity'); ylabel('x');
ax = gca; ax.FontSize = 14; axis tight;

% Just some of the range 
figure(2);
plot(fb(150:300), x(150:300),'m');
title('Optimal stopping boundary for each node x');
xlabel('Time to maturity'); ylabel('x');
ax = gca; ax.FontSize = 14; axis tight;
    



