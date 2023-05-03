%%
%% Process for in-vitro images
%% 
%% Extract z-stacked image from the .lif files obtained from the microscope

% There are differen way to extract the z stacked image, just be sure to
% obtain a 2D image that projects the z stacked obtained.

%% Once saved the stacked image, we manually segment each organoid
% Note: you may find other segmentation methods that can work for some
% organoids. No matter how you segment the organoids you can still use the
% code to detect/count crypt-like structures.

load('stackedim_serie_02.mat')
% I = imread('C:\Users\sm16476\Documents\SimpleCryptCount_Project-DEMO\Extract_z_stack\Organoid Day 3 .lif_folder\stackedim_serie_02.png')

I = im_stack;
imshow(I)
imageSegmenter(I)

%% After exporting the mask ('BW') and the masked image ('maskedImage')
binaryImage = BW;

save('Org2_example', 'binaryImage')

%% We can use that BW on the CountingCrypts_wCircularityFun.m

%Crypt parameters:
Input_min_area = 0.0666;
Input_max_area = 0.2736;
Input_min_arcLength = 0.1466;

[NumCrypts Circularity] = CountingCrypts_wCircularityFun ('In vitro', 'Org1_example', 7.8, [Input_min_area, Input_max_area, Input_min_arcLength])

%% Adjust the parameters for your own data

% Set an initial value for all parameters
Input_min_area = 0.0666;
Input_max_area = 0.2736;
Input_min_arcLength = 0.1466;
fourier_harmonic_term = 7;

x = [Input_min_area, Input_max_area, Input_min_arcLength, fourier_harmonic_term];

error_D3 = simple_objective_manualSeg_D3(x)
