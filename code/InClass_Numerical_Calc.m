clc
clear 
close all

%Initial condition
N0 = 1; %number of bacteria at t = 0
k = 1/30; % every 30 min
dt = 15; % shoulb be smaller than 30 mins 


N(1) = N0; %initial value of N
time(1) = 0; % initial value of time

for i = 2: 50
    
    N(i) = N(i-1) + k*dt*N(i-1);
    time(i) = (i-1)*dt;
    
end

figure(1)
semilogy(time,N)
xlabel('time (min)')
ylabel('Number of cells')

% hold on
% 
% k = 1/50; % every 50 min bacteria 
% 
% for i = 2: 50
%     
%     N(i) = N(i-1) + k*dt*N(i-1);
%     time(i) = (i-1)*dt;
%     
% end
% 
% semilogy(time,N,'r')
% 
% 
% k = 1/15; % every 15 min bacteria 
% 
% for i = 2: 50
%     

%% upload the data

data = csvread('Growth_area.csv');
time = data(1,:);
area = data(2,:);

figure(2)
plot(time,area,'.k')

xlabel('time (min)')
ylabel('area (a.u.)')

figure(3)
semilogy(time,area,'.k')

xlabel('time (min)')
ylabel('area (a.u.)')

%% Ki squared mathod

% 0.0001 = 1e-4

kRange  = 0.01:1e-4:0.5;


N0 = area(1);

for i = 1:length(kRange)
    
   error(i) = sum((log(area) - (log(N0) + kRange(i)*time)).^2);
    
end

figure(4)
semilogy(kRange,error,'o')
xlabel('k(1/min)')
ylabel('error')


%% find the best fit

kFit = 0.0264;

t_doubling = log(2)/kFit; % t_doubling = 26.2556 min

figure(5)
semilogy(time, area,'ro')
hold on

best_fit = N0*exp(kFit*time);
semilogy(time,best_fit)
hold off

xlabel('time (min)')
ylabel('Cell area')

legend('Experimental area','Area from the best fit model')









