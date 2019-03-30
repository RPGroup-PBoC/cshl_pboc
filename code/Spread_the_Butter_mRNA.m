clc 
clear 
close all

%% Model parameters
r = 3;            % 1/min rate of transcription of mRNA
gamma = 1/3;      % 1/min rate of decay of mRNA

Nmax = 40; % max number of mRNA
m = 0:1:Nmax; % number of mRNA

dt = 0.01; % min
TotalTime = 50; %min

time = 0:dt:TotalTime; % time vector

%% Initialize the prob
p = zeros(length(m),length(time));

% initial condition
p(1,1) = 1;% probability at t = 0

for t = 2:length(time)
    
    % zero number of mRNA
    
    p(1,t) = p(1,t-1) - r*dt*p(1,t-1) + gamma*dt*(2-1)*p(2,t-1);
    
    % max number of mRNA
    
    p(length(m),t) = p(length(m),t-1) + r * dt * p(length(m)-1,t-1)...
        -gamma * dt * Nmax * p(length(m),t-1);
    
    for n = 2:length(m)-1
        
        p(n,t) = p(n,t-1) + r*dt*p(n-1,t-1) - r*dt*p(n,t-1) ...
            + gamma*dt*n*p(n+1,t-1) - gamma*dt*(n-1)*p(n,t-1);
        
    end
       

end

figure(1)

for i = 1:length(time)
    bar(m,p(:,i))
    ylim([0 1])
    xlabel('mRNA number')
    ylabel('Probability')

    drawnow
    pause(0.1)  
end

figure(2)

bar(m,p(:,length(time)))

hold on

plot(m,(r/gamma).^m * exp(-r/gamma)./factorial(m),'.-r')

hold off

xlabel('mRNA number')
ylabel('Probability')

legend('Simulation','Poisson')