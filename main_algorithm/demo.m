addpath('../DAISY/');   % DAISY feature extraction
addpath('../SLIC-Superpixel-master/'); % SLIC superpixel algorithm

%% Loading RF Model
warning('off', 'MATLAB:class:mustReturnObject');
load('model/Mdl_25_spold2.mat');

foreground_path ='..\Dataset_test\ISTD_small_A\';
mask_path = '..\Dataset_test\ISTD_small_mask\';
result_path = 'results_out\'; % Change to your preferred path

image_search   = '*.png';
foregrounds   = dir ([foreground_path, image_search]);
masks           = dir ([mask_path,image_search]); 
write_img     = true;
if ~exist(result_path, 'dir')
   mkdir(result_path)
end

%% Parameters
params.plot_figures = 0;  % whether to plot results (0 or 1)
params.model = Mdl; % RF Model path
params.verbose = 0;   %Verbosity: 0 (no print)/ 1(print message)/ 2(print more messages)

%% START OF THE ALGORITHM
for index = 1:size(foregrounds,1) 
    t_start = tic();
    foreground_filename = foregrounds(index).name;
    mask_filename = masks(index).name;    
    
    image_name = strsplit(foreground_filename, '.');
    image_name = char(image_name(1));
    fprintf('========================START========================\n'); 
    fprintf('Processing image %s (%d/%d)\n',image_name,index,size(foregrounds,1)); 
    
    % create filename
    write_path = [result_path,image_name,'_output'];
    
    % load images
    foreground = im2double(imread([foreground_path,foreground_filename]));
    mask = im2double(imread([mask_path,mask_filename]));   
    mask = double(mask(:,:,1)>0.5);
    [ shadow_free  ] = shadow_removal_icme22( foreground, mask, params);
       
    t_end = toc(t_start);
    if(write_img)
        imwrite(shadow_free, [write_path, '.png']);
    end
    
    fprintf('Elapsed time: %2.5f secs\n', t_end);
end