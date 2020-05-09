% Reads in data
format long;
value = csvread('EllipticIneqPlot.csv');

nodes = value(1,1:end-1);  % x values
approx = value(2, 1:end-1); % u approximation
exact = value(3, 1:end-1); % exact u values

% Create plot
figure(1);
plot(nodes,approx,'m',nodes,exact, 'c', 'LineWidth',2);
xlabel('x'); ylabel('u(x)'); legend('Approximation','Unconstrained');
title('Solution to elliptic inequality PDE problem'); 
ax = gca; ax.FontSize = 14; axis tight;
