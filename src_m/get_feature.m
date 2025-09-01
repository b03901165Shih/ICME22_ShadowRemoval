% 15+15-3 features
function [X] = get_feature(Lin_feat, Lout_feat, cov_feat, num_samp_Lin, u_Lout, offset)

    X = zeros(num_samp_Lin, max(u_Lout(:)), 14);
    X_sort = zeros(num_samp_Lin, max(u_Lout(:)), size(X,3)-3);

    % Find the feature vector for every pair
    for i = 1:num_samp_Lin
        % Positive match
        for j = 1:max(u_Lout(:))  
            % average of rgb
            X(i,j, 1:3) = Lin_feat{i,1};
            % rgb ratios
            X(i, j, 4:6) = (Lin_feat{i,1}+offset)./(Lout_feat{j,1}+offset);
            % earth mover distance of luminance hist
            X(i, j, 7) = pdist2_mod(Lin_feat{i,2}, Lout_feat{j,2}, 'emd');
            % median a, b, offset
            X(i, j, 8:9) = [Lin_feat{i,6}- Lout_feat{j,6}, Lin_feat{i,7}- Lout_feat{j,7}];

            %grad mag distance  
            X(i, j, 10) = pdist2_mod(Lin_feat{i,4}, Lout_feat{j,4}, 'emd');
            %grad angle distance  
            X(i, j, 11) = pdist2_mod(Lin_feat{i,5}, Lout_feat{j,5}, 'emd');
            %daisy distance         
            X(i, j, 12) = sqrt(sum((Lin_feat{i,3}-Lout_feat{j,3}).^2));
            %gabor distance        
            X(i, j, 13) = sqrt(sum((Lin_feat{i,8}-Lout_feat{j,8}).^2));
            %spatial distance        
            %X(i, j, 14) = sqrt(sum((Lin_feat{i,9}-Lout_feat{j,9}).^2));
            
            %covariance mapped | average of rgb        
            %X(i, j, 15:17) = cov_feat{i,j,1};
            %covariance mapped | correlation   
            X(i, j, 14) = cov_feat{i,j,2};
            %covariance mapped | earth mover distance of luminance hist   
            %X(i, j, 19) = pdist2_mod(cov_feat{i,j,3}, Lout_feat{j,2}, 'emd');
        end       
    
        %Sort:
        
        [~, X_sort_idx]=sort(squeeze(X(i,:,4:end)));
        [~, X_sort_idx]=sort(X_sort_idx); % inverse indexing
        X_sort(i,:,:) = X_sort_idx/(max(u_Lout(:)) );
        
        
        %Normalized feature
        %X_tmp = squeeze(X(i,:,4:end));
        %X_sort(i,:,:) = bsxfun(@rdivide, bsxfun(@minus,X_tmp, mean(X_tmp,1)), std(X_tmp,1));
    end
    X = cat(3, X, X_sort);
end