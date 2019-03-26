
clc % to clean the display on command window
close all % to close all the MATLAB windows 
clear % to clear all the variables in MATLAB before you start your calculations


a = [4 5 6];
b = [1;2;3];
c = [1 2 3; 4 5 6];

a(1);

c>3; %binary matix , logical matix, only has values of 1 and zero

c(c>3); 

%% Image Analysis

cd colony_growth % change directory to the folder of images
image = imread('EcoliGrowth14.tif'); % matlab loads the image in a matrix

figure(1)
imshow(image) % matlab shows the image that it read



imtool(image); % to inspect the value of pixels in the image


threshold = 100;

image_threshold = image < threshold;

figure(2)
imshow(image_threshold,[])

area = sum(image_threshold(:));% sum all the values in all rows and all clms %sum(sum(image_threshold))
%area_2 = sum(sum(image_threshold));


%% Load the images 

images = dir('*.tif');
threshold = 100;

for i = 1:length(images)
    image = imread(images(i).name);
    image_threshold = image<threshold;
    area = sum(sum(image_threshold));
    area_list(i) = area;
end

figure(3)
plot(area_list)
xlabel('number of images')
ylabel('area for each image')

time_btw_frame = 5; %min
time = 0:time_btw_frame:(length(images)-1)*time_btw_frame;

figure(4)
plot(time,area_list)
xlabel('time (min)')
ylabel('area (a.u.)')


figure(5)
semilogy(time, area_list)
xlabel('time (min)')
ylabel('area (a.u.)')

figure(6)
plot(time, log(area_list))
xlabel('time (min)')
ylabel('log(area) (a.u.)')

