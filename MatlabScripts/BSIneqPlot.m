% Reads in data
format long;
value = csvread('BSIneqPlot.csv');

nodes = value(1,1:end-1); 
initU = value(2, 1:end-1);
uT2 = value(3, 1:end-1);
uT2Exact = value(4, 1:end-1);
uT = value(5, 1:end-1);
uTExact = value(6, 1:end-1);

% Create plots
figure(1);
plot(nodes, initU, 'm');
title('American put option value at t=0 (option payoff)');
xlabel('x'); ylabel('u(x)');
ax = gca; ax.FontSize = 14; axis tight;

figure(2);
plot(nodes,uT2,'m',nodes,uT2Exact, 'c', 'LineWidth',2);
xlabel('x'); ylabel('u(x)'); legend('Approximation','European');
title('American put option value at t=T/2'); 

figure(3);
plot(nodes,uT,'m',nodes,uTExact, 'c', 'LineWidth',2);
xlabel('x'); ylabel('u(x)'); legend('Approximation','European');
title('American put option value at t=T (Option price)'); 
