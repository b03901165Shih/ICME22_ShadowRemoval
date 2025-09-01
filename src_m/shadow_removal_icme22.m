function [ shadow_free  ] = shadow_removal_icme22( foreground, mask, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if(~isfield(params, 'num_iters'))
        params.num_iters = 1;
    end
    if(~isfield(params, 'optimization'))
        params.optimization = 1;
    end
    if(~isfield(params, 'aggregate_cluster'))
        params.aggregate_cluster = 1;
    end    
    if(~isfield(params, 'boundary_opt'))
        params.boundary_opt = 1;
    end


    mask_dilate = imdilate(mask, strel('disk', 4));
    
    %% Main algo start
     [candidates,  L_in, L_out ] = find_cluster(params, foreground, mask_dilate );
    
    %% Clusters Aggregation
    [ clustered_L_in, clustered_L_out ] = cluster_agg( params, foreground, L_in, L_out);
    
    %% Predict positive/negative matching using pre-trained model
    if(params.verbose>0)
        fprintf('Testing feature predictions\n');
    end
    [ Y ] = get_testing_pred( foreground, clustered_L_in, clustered_L_out, params.model );
    
    %% Finding mapping functions 
    if(params.verbose>0)
        fprintf('Recovering mapping functions using bayesian\n');
    end
    [ remapped_Im ] = recover_mapping_from_prediction( foreground, clustered_L_in, clustered_L_out, Y );
    remapped_f = remapped_Im.*repmat(mask_dilate==1,[1,1,3])+foreground.*repmat(mask_dilate==0,[1,1,3]);    
    
    if(params.plot_figures)
        figure, imshow(remapped_f);
    end
    
    %% Optimize mapping functions     
    if(params.optimization)
        if(params.verbose>0)
            fprintf('Minimize energy by optimizing mapping functions\n');    
        end
        [ f_t, ~, E(1) ] = func_energy_min(params, foreground, remapped_f, clustered_L_in, clustered_L_out, Y );
        [ remapped_Im ] = recover_mapping_from_prediction( foreground, clustered_L_in, clustered_L_out, Y , f_t);

        remapped_f = remapped_Im.*repmat(mask_dilate==1,[1,1,3])+foreground.*repmat(mask_dilate==0,[1,1,3]);
        candidates(:,:,:,1) = remapped_f;

        if(params.plot_figures)
            figure, imshow(remapped_f);
        end

        %% Post-processing (iterative optimize) 
        for i = 1:params.num_iters
            foreground_trans = remapped_f;%(foreground+remapped_f)*0.5;
            if(params.plot_figures)
                figure, imshow(foreground_trans);
            end
            [ Y ] = get_testing_pred( foreground_trans, clustered_L_in, clustered_L_out, params.model );

            if(params.verbose>0)
                fprintf('Recovering mapping function using bayesian\n');
            end
            [ remapped_Im ] = recover_mapping_from_prediction( foreground_trans, clustered_L_in, clustered_L_out, Y );
            remapped_f = remapped_Im.*repmat(mask_dilate==1,[1,1,3])+foreground_trans.*repmat(mask_dilate==0,[1,1,3]);        
            if(params.plot_figures)
                figure, imshow(remapped_f);
            end

            if(params.verbose>0)
                fprintf('Minimize energy by optimizing mapping functions\n');    
            end
            [ f_t, ~, E(i+1) ] = func_energy_min(params, foreground_trans, remapped_f, clustered_L_in, clustered_L_out, Y );
            [ remapped_Im ] = recover_mapping_from_prediction( foreground_trans, clustered_L_in, clustered_L_out, Y , f_t);     
            remapped_f = remapped_Im.*repmat(mask_dilate==1,[1,1,3])+foreground_trans.*repmat(mask_dilate==0,[1,1,3]);    
            if(params.plot_figures)
                figure, imshow(remapped_f);
            end
            candidates(:,:,:,i+1) = remapped_f;
        end

        %[~, idx] = min(E);
        shadow_free = squeeze(candidates(:,:,:,end));
    end
    
    %% Boundary inpainting
    if(params.boundary_opt)
        se_bw = strel('disk',4);
        mask_dilate = imdilate(mask, strel('disk', 2));
        BW = (imdilate(mask_dilate,se_bw) > mask_dilate) | (imerode(mask_dilate,se_bw) < mask_dilate);
        %BW = (imdilate(mask,se_bw) > mask) | (imerode(mask,se_bw) < mask);
        shadow_free = inpaintExemplar(shadow_free, BW);
        if(params.plot_figures)
            figure, imshow(shadow_free);
        end
    end
    
end