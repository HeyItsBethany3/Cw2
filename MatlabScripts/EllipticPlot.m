% Reads in data
format long;
value = csvread('EllipticPlot.csv');

nodes = value(1,1:end-1); 
approx = value(2, 1:end-1);
exact = value(3, 1:end-1);

% Create plot
figure(1);
plot(nodes,approx,'m',nodes,exact, 'c', 'LineWidth',2);
xlabel('x'); ylabel('u(x)'); legend('Approximation','Exact');
title('Solutiion to elliptic PDE problem with n=8'); % change this
ax = gca; ax.FontSize = 14; axis tight;


