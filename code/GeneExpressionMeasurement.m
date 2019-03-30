%Measure the fluorescence per cell as a means to obtain
%the level of gene expression and eventually calculate
%the fold-change in gene expression.

%Approach:
%1) Find the cells in the phase contrast image by
%   thresholding.
%2) Figure out which pixels in the fluorescence image
%   correspond to each cell.
%3) Obtain the fluorescence per cell.

%Load the phase contrast image
ImPhase=imread('PhaseConstrastImage.tif');
%Display the phase contrast image
imshow(ImPhase)
%Matlab doesn't know how to display a 16-bit image which has (65k)
%levels of grey. It only knows how to display 256 levels of grey.
%One way around this is to ask Matlab to assign 0 to the dimmest
%pixel in the image and 255 to the brightest pixel.
imshow(ImPhase,[])

%Segment the cells by thresholding. First, I need to find a threshold.
%imtool(ImPhase,[])  %I'm commenting this out so that it doesn't
                     %load a new imtool image every time I run the code
Threshold=3500;

%Take threshold
ImThresh=ImPhase<Threshold;
imshow(ImThresh)

%We did a relatively good job at segmenting the cells, but we also
%captured some speckles / crap. We want to get rid of this by
%taking into account the fact that cells have a given area.
%First, we need to assign each cell a unique identity using bwlabel.
ImLabel=bwlabel(ImThresh);
imshow(ImLabel,[])
%Calculate the area of each region
ImProps=regionprops(ImLabel,'area');
%ImProps is a structure array containing one field for the Area.
%If I want the area corresponding to the 15th region, I do
ImProps(15).Area
%Let's get a vector with all the areas
Areas=[ImProps.Area];
%Plot the histogram of the areas using 50 bins
histogram(Areas,50)
xlabel('area (pixels)')
ylabel('number of regions')

%We are going to filter out any regions with an area smaller
%than 50 pixels:
%1) Go through each region and ask whether its area is higher
%   than the treshold.
%2) If yes, I need to isolate that region and add it to
%   blank image.

%How do you isolate a specific region in our labelled image?
%Let's show, for example, only region number 10
imshow(ImLabel==10)

%Define the area threshold
AreaThresh=50;
%Define a blank image
ImNew=zeros(size(ImLabel));
%Loop through all regions
figure(1)
for i=1:length(Areas)
   %Is the area larger than AreaThresh?
   if Areas(i)>AreaThresh
       ImNew=ImNew+(ImLabel==i);
       %imshow(ImNew)
       %pause(0.5)
   end
end

%To calculate the fluorescence, I need to multiply a mask
%corresponding to each cell with the fluorescence image.

%Load the fluorescence image
ImFluo=imread('FluorescenceImage.tif');

%I need to relabel my image ImNew
ImLabel2=bwlabel(ImNew);

%Isolate the fluorescence of cell number 10
imshow(immultiply(ImLabel2==10,ImFluo),[])
%Calculate the fluorescence of this cell
sum(sum(immultiply(ImLabel2==10,ImFluo)))

%Now do it for all cells
for i=1:max(max(ImLabel2))
   %Generate a mask from the thresholded image
   ImMask=(ImLabel2==i);
   %Isolate the fluorescent pixels corresponding to the mask
   FluoMask=immultiply(ImMask,ImFluo);
   %Calculate and save the fluorescence
   CellFluo(i)=sum(sum(FluoMask));
end

histogram(CellFluo,20)
xlabel('cell fluorescence')
ylabel('number of cells')












