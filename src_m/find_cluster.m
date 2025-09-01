function [candidates,  clustered_L_in, clustered_L_out ] = find_cluster(params, foreground, mask )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    candidates = zeros(size(foreground, 1), size(foreground, 2), size(foreground, 3), params.num_iters+1);
        
    if(params.plot_figures)
        figure, imshow(foreground);
    end
    %mask_dilate = mask;
    
    edge_width = 100;
    mask_out = imdilate(mask, strel('disk',edge_width))-mask;
    mask_in  = mask-imerode(mask, strel('disk',edge_width));
    
    %% Clustering
    if(params.verbose>0)
        fprintf('Extracting clusters from masked image\n');
    end
    [ clustered_L_in, ~ ] = superpixel_clustering( foreground, mask_in, 1024 );
    [ clustered_L_out, ~ ] = superpixel_clustering( foreground, mask_out, 1024 );        

end

function [ L_in, L_out ] = superpixel_clustering( foreground, mask_dilate, K )

    num_suppix = K;
    L  = mexGenerateSuperPixel(foreground, num_suppix, 10, mask_dilate)+1;

    %basicC = basic;
    %basicL = rgb2hsv(basic);
    %basicL = basicL(:,:,1);

    %sigma_c = std(reshape(basicC,[],3)); sigma_l = std(basicL(:));

    %% reodering Lin and L_out
    L_in = L.*mask_dilate;
    L_in_ele = unique(L_in);
    tmpL_in = L_in;
    for i = 1:size(L_in_ele,1)
        tmpL_in(L_in==L_in_ele(i)) = i; 
    end
    L_in = tmpL_in-1;

    L_out = L.*(1-mask_dilate);
    L_out_ele = unique(L_out);
    tmpL_out = L_out;
    for i = 1:size(L_out_ele,1)
        tmpL_out(L_out==L_out_ele(i)) = i; 
    end
    L_out = tmpL_out-1;

end
