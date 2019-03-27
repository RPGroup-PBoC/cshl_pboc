clc
close all
clear

%rand() % uniform distribution btw (0,1)

% number of walker 1

position = 0; % position at t = 0
step_size = 1;
prob_right = 0.5;
prob_left = 1 - prob_right;

coin_flip = rand();

if coin_flip <= prob_left
    position = position - step_size;
else
    position = position + step_size;
end

n_steps = 50;


%% for one walker
position(1) = 0;

for i = 2:n_steps
    
    coin_flip = rand();
    if coin_flip <= prob_left
        position(i) = position(i-1) - step_size;
    else
        position(i) = position(i-1) + step_size;
    end
    
end

time = 1:1:n_steps;

figure(1)
plot(time,position)
xlabel('time')
ylabel('position')
title('1 walker')


%% multiple walkers

n_steps = 50; 
n_walkers = 4000;
step_size = 1;
prob_left = 0.5;
prob_right = 1- prob_left;

position = zeros(n_walkers,n_steps); % to intialize the position matix 

for i = 1:n_walkers
    for j = 2:n_steps
        coin_flip = rand();
        if coin_flip<= prob_left
            position(i,j) = position(i,j-1) - step_size;
        else
            
            position(i,j) = position(i,j-1) + step_size;
        end
        
    end
end


time = 1:1:n_steps;

% figure(2)
% 
% for i = 1:n_walkers
%     plot(time,position(i,:)) % plot each walker for all the n steps
%     ylim([-20 20]) %range of y axis 
%     xlabel('time')
%     ylabel('position')
%     title('20 walker')
%     hold on
%     pause(0.3)
% end
% 
% hold off


%% Calculate MSD

disp = zeros(n_steps,n_walkers);

for i = 2:n_steps  
    
    for j = 1:n_walkers
        disp(i,j) = position(j,i) - position(j,1);
    end
   
end

disp_sq = disp.^2;

mean_disp = zeros(1,n_steps);

for i = 1:n_steps
    
    mean_disp(i) = mean(disp_sq(i,:)); % taking the mean of all the walkers in step i
    
end


figure(3)
time = 0:1:n_steps-1;
plot(time,mean_disp)

hold on
expected_mean_disp = time;
plot(time,expected_mean_disp,'r')

xlabel('time')
ylabel('MSD')


hold off