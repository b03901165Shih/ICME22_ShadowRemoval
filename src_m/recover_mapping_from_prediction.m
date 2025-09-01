function [ remapped_Im ] = recover_mapping_from_prediction( foreground, clustered_L_in, clustered_L_out, Y, f_t )

match_idx = zeros(size(Y));
match_score = zeros(size(Y));
for i = 1:size(Y,1)
    [match_score(i,:), match_idx(i,:)]=sort(Y(i,:), 'descend');
end
se = ones(3);

%%
    count = 1;
    %% Find Aj, Bj, muj, sigmaj for all clusters (stored in models)
    for i = 1:size(Y, 1)
        inner_mask = imclose(clustered_L_in==i, se);
        RGB_shadow = reshape(foreground(repmat(inner_mask , [1,1,3])), [], 3);
        if nargin < 5 || isempty(f_t)
            threshold = min(0.5,match_score(i, 1)); % match_score(i, 1); %
            if(match_score(i,1)<threshold), continue; end
            outer_mask = imclose(ismember(clustered_L_out,match_idx(i,match_score(i, :)>=threshold)), se);%match_score(i,:)>=0.5
            RGB_lit = reshape(foreground(repmat(outer_mask, [1,1,3])), [], 3);
            [ model.A, model.B ] = find_covariance_mapping( RGB_lit, RGB_shadow); %find_covariance_mapping
        else
            model.A = f_t{i,1};
            model.B = f_t{i,2};
        end
        model.mu = mean(RGB_shadow, 1);
         c_RGB_shadow = bsxfun(@minus, RGB_shadow, model.mu);
         model.sigma = (((c_RGB_shadow)'*(c_RGB_shadow))/size(RGB_shadow,1))+eps * eye(3);
         model.num = size(RGB_shadow,1);
         models(count) = model;
         count = count +1;
    end
    
    %% Also force lit pixel to be identitcal
    %%{
    for i = 1:size(Y, 2)
        outer_mask = imclose(clustered_L_out==i, se);
        RGB_lit = reshape(foreground(repmat(outer_mask, [1,1,3])), [], 3);
        model.mu = mean(RGB_lit, 1);
        c_RGB_lit = bsxfun(@minus, RGB_lit, model.mu);
        model.A = eye(3);
        model.B = [0,0,0];
        model.sigma = (((c_RGB_lit)'*(c_RGB_lit))/size(RGB_lit,1))+eps * eye(3);
        model.num = size(RGB_lit,1)*0.01;
         models(count) = model;
         count = count +1;
    end
    %%}
    
    %% Construct the mapped image using the mapping function of each cluster
    weight_map = zeros(length(models), size(foreground,1), size(foreground,2));
    for cluster_id = 1: size(models, 2)
        weight_map(cluster_id,:,:) = models(cluster_id).num* ...
                reshape(mvnpdf(reshape(foreground,[],3),models(cluster_id).mu,models(cluster_id).sigma), size(foreground,1), size(foreground,2))+eps;
    end
    weight_sum = sum(weight_map, 1);
    for cluster_id = 1: size(models, 2)
        weight_map(cluster_id,:,:) = weight_map(cluster_id,:,:)./weight_sum;
    end    
    remapped_Im = zeros(size(foreground));
    for cluster_id = 1: size(models, 2)
        A = models(cluster_id).A;
        B = models(cluster_id).B;
        w = repmat(squeeze(weight_map(cluster_id,:,:)), 1, 1, 3);
        remapped_Im = remapped_Im+w.*(reshape(reshape(foreground, [], 3)*A+repmat(B, size(foreground, 1)*size(foreground, 2), 1), size(foreground))); 
    end
    %figure, imshow(remapped_Im)



end



%%
% Mapping the covariance and mean of target_vec to that of ref_vec
% ref_vec, target_vec: n x p matrix, n: samples, p: dimension
% A: p x p matrix
% B: 1 x p
% remapped_target_vec = (target_vec*A+repmat(B, N2, 1))

%{
function [ A, B] = find_covariance_mapping( ref_vec, target_vec)

% numbers of data
N1 = size(ref_vec, 1);
N2 = size(target_vec, 1);

%Matching mean only
mean_ref = mean(ref_vec, 1);
mean_target = mean(target_vec, 1);

% Cholesky decomposing the covariance matrix to upper triangle matrix 
coeff_ref_vec = (ref_vec-repmat(mean_ref, N1, 1))'*(ref_vec-repmat(mean_ref, N1, 1))/(N1-1)+1e-10.*eye(3);
coeff_target_vec = (target_vec-repmat(mean_target, N2, 1))'*(target_vec-repmat(mean_target, N2, 1))/(N2-1)+1e-10.*eye(3);

U_ref = chol(coeff_ref_vec); 
U_target = chol(coeff_target_vec);
A = real(pinv(U_target)*U_ref);

%[L_ref, D_ref] = ldl(coeff_ref_vec);
%[L_target, D_target] = ldl(coeff_target_vec);
%A = real(pinv(sqrt(D_target)'*L_target')*(sqrt(D_ref)'*L_ref'));

B = mean_ref-mean_target*A;
end
%}
