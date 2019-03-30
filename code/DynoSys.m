% for the mutual repression genetic switch, we want to plot the phase
% portrait. the means that we want to calculate the velocity vector for
% each possible concentration value of R1 and R2.
clc
close all
clear 
% dR1/dt = r( 1/(1+(R1/Kd)^2) 0 - gamma * R1
% dR2/dt = r( 1/(1+(R2/Kd)^2) ) - gamma * R2
% SO WE HAVE:
% R2 = r/gamma * 1/(1+(R1/Kd)^2)
%% model parameters
r = 20;         % rate of protein production (protein/min)
gamma = 1/30;   % rate of degradation in 1/(min)
Kd = 200;       % dissociation constant in units of number of proteins.


% first, plot two node points.
Range = 600; %max amount of protein 
R = linspace(0, Range,65); % divide interval btw 0 and Range to 65 smaller intervals


figure(1); 
% Nullcline R2
p(1) = plot(R, r/gamma*1./(1+(R./Kd).^2), 'r', 'linewidth', 1.4);
hold on;
p(2) = plot(r/gamma*1./(1+(R./Kd).^2), R, 'color', [0,0.6,0], 'linewidth', 1.4);
xlabel('R1', 'fontsize', 20);
ylabel('R2', 'fontsize', 20);


% to plot the arrows on the phase portrat, we need to create a matrix tht
% has all possible combinations of R1 and R2 values. this is done using
% meshgrid.
[R1M, R2M] = meshgrid(R(1:2:end), R(1:2:end));

dR1M = -gamma.*R1M+r./(1+(R2M/Kd).^2); %ODE number 1
dR2M = -gamma.*R2M+r./(1+(R1M/Kd).^2); %ODE number 2
hold on; 
p(3) = quiver(R1M, R2M, dR1M, dR2M, 1.5, 'linewidth',1);

h=legend(p(1:2),'dR2/dt = 0','dR1/dt = 0'); set(h,'fontsize',15)
axis equal 
xlim([1 Range])
set(gca,'FontSize',18)
% now we want to simulate the dynamics of R1 and R2 given some initial
% condition.
% step1, define initial R1 and R2;
N = 10; % number of random initial conditions
R1Ini = rand(N,1).*Range;
R2Ini = rand(N,1).*Range;
% then loop over many many time steps...
R1 = R1Ini;
R2 = R2Ini;
for t = 1:500 % number of steps to follow the trajectory
    % at each tiem point, calculate the local evolution direction
    dR1 = -gamma.*R1+r./(1+(R2/Kd).^2);
    dR2 = -gamma.*R2+r./(1+(R1/Kd).^2);
    R1 = R1 + dR1; % integration is at dt = 1 
    R2 = R2 + dR2;
    %figure(1); 
    plot(R1, R2, 'm.');%,'Visible','off','displayname','data'); 
    drawnow
    pause(0.2)
    legend off
end



