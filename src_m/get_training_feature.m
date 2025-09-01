function [ X, Y ] = get_training_feature( foreground, chi_match, chi_distance, clustered_L_in, clustered_L_out  )

u_Lin = unique(clustered_L_in);
u_Lout = unique(clustered_L_out);

%RGB ratio map
offset = 1e-3;
foreground_rgb = cat(3, mean(foreground, 3), ((foreground(:,:,1)+offset)./(foreground(:,:,3)+offset)), ((foreground(:,:,2)+offset)./(foreground(:,:,3)+offset)));

%%{
gabor_f = OwnGetGaborResponse(rgb2gray(foreground));
gabor_f = bsxfun(@minus, gabor_f, min(gabor_f,[],3));
gabor_f = bsxfun(@rdivide, gabor_f, max(gabor_f,[],3));
%%}

[angles, magnit] = gradient_profile(mean(foreground,3));
magnit = magnit/max(magnit(:));
angles = (angles+pi)/2/pi;

pad_r = 15; pca_dim = 25;
dzy = compute_daisy(padarray(foreground,[pad_r, pad_r], 'replicate'));
[~, score] = pca(dzy.descs);
DAISY_f = permute(reshape(score(:,1:pca_dim), size(foreground,2)+2*pad_r, size(foreground, 1)+2*pad_r,[]),[2,1,3]);
DAISY_f = DAISY_f(pad_r+1:end-pad_r,pad_r+1:end-pad_r,:);
DAISY_f = bsxfun(@minus, DAISY_f, min(DAISY_f,[],3));
DAISY_f = bsxfun(@rdivide, DAISY_f, max(DAISY_f,[],3));

% Luminance 
I_lab = rgb2lab(foreground);
I_lab(:,:,1) = I_lab(:,:,1)/100;
I_lab(:,:,2) = (I_lab(:,:,2)+100)/(200); 
I_lab(:,:,3) = (I_lab(:,:,3)+100)/(200); 
edges = 0:0.02:1;
l_edges = 0:0.1:1;

%Distance
[Y_coord, X_coord] = meshgrid(1/size(foreground, 2):1/size(foreground, 2):1, 1/size(foreground,1):1/size(foreground,1):1);


%% Extract covariance features from regions
% Covariance Mapped result
RGB_shadow = cell(max(u_Lin),1); RGB_lit = cell(max(u_Lout),1);
for i = 1:max(u_Lin(:)), RGB_shadow(i, 1) = {reshape(foreground(repmat(clustered_L_in==i,[1,1,3])),[],3)}; end
for i = 1:max(u_Lout(:)),  RGB_lit(i, 1) = {reshape(foreground(repmat(clustered_L_out==i,[1,1,3])),[],3)}; end

cov_feat = cell(max(u_Lin), max(u_Lout), 3);
for i = 1:max(u_Lin(:))
    for j = 1:max(u_Lout(:))
        [ A, B ] = find_covariance_mapping( RGB_lit{j}, RGB_shadow{i}); %find_covariance_mapping
        RGB_relit = bsxfun(@plus, RGB_shadow{i}*A, B);
        RGB_relit_rgb = cat(3, mean(RGB_relit, 2), ((RGB_relit(:,1)+offset)./(RGB_relit(:,3)+offset)), ((RGB_relit(:,2)+offset)./(RGB_relit(:,3)+offset)));
        cov_feat(i,j, 1) = {mean(RGB_relit_rgb, 1)};%average rgb 
        cov_feat(i,j, 2) = {corr2(RGB_relit, RGB_shadow{i})}; % correlation to original pixels
        cov_feat(i,j, 3) = {histcounts(RGB_relit_rgb(:,1), edges,'Normalization', 'probability')}; % Histogram of mean intensity
    end
end

%% Extract features from regions
num_samp_Lin = max(u_Lin(:));
Lin_feat = cell(num_samp_Lin, 9);
Lout_feat = cell(max(u_Lout(:)), 9);

for i = 1:num_samp_Lin    
    L_in_mask = (clustered_L_in==i);
    lab_extract = reshape(I_lab(repmat(L_in_mask, [1,1,3])),[],3);
    daisy_extract = reshape(DAISY_f(repmat(L_in_mask, [1,1,size(DAISY_f, 3)])),[],size(DAISY_f, 3));
    gabor_extract = reshape(gabor_f(repmat(L_in_mask, [1,1,size(gabor_f, 3)])),[],size(gabor_f, 3));
    Lin_feat(i,1) = {mean(reshape(foreground_rgb(repmat(L_in_mask, [1,1,3])),[],3), 1)}; %average rgb
    Lin_feat(i,2) = {histcounts(lab_extract(:,1), edges,'Normalization', 'probability')}; %luminance histogram
    
    Lin_feat(i,3) = {mean(daisy_extract)}; %average daisy
    Lin_feat(i,4) = {histcounts(magnit(L_in_mask), l_edges,'Normalization', 'probability')}; % maginutude histogram
    Lin_feat(i,5) = {histcounts_w(angles(L_in_mask), magnit(L_in_mask), l_edges)}; % maginutude histogram
    Lin_feat(i,6) = {median(lab_extract(:,2))}; % average a
    Lin_feat(i,7) = {median(lab_extract(:,3))}; % average b
    Lin_feat(i,8) = {mean(gabor_extract)}; %average gabor
    Lin_feat(i,9) = {[mean(X_coord(L_in_mask)), mean(Y_coord(L_in_mask))]}; %average distance
end
for i = 1:max(u_Lout(:))
    %if(~ismember(i,unique_out)), continue; end
    L_out_mask = (clustered_L_out==i);
    lab_extract = reshape(I_lab(repmat(L_out_mask, [1,1,3])),[],3);
    daisy_extract = reshape(DAISY_f(repmat(L_out_mask, [1,1,size(DAISY_f, 3)])),[],size(DAISY_f, 3));
    gabor_extract = reshape(gabor_f(repmat(L_out_mask, [1,1,size(gabor_f, 3)])),[],size(gabor_f, 3));
    Lout_feat(i,1) = {mean(reshape(foreground_rgb(repmat(L_out_mask, [1,1,3])),[],3), 1)}; %average rgb
    Lout_feat(i,2) = {histcounts(lab_extract(:,1), edges,'Normalization', 'probability')}; %luminance histogram
    
    Lout_feat(i,3) = {mean(daisy_extract)}; %average daisy
    Lout_feat(i,4) = {histcounts(magnit(L_out_mask), l_edges,'Normalization', 'probability')}; % maginutude histogram
    Lout_feat(i,5) = {histcounts_w(angles(L_out_mask), magnit(L_out_mask), l_edges)}; % angle histogram
    Lout_feat(i,6) = {median(lab_extract(:,2))}; % average a
    Lout_feat(i,7) = {median(lab_extract(:,3))}; % average b
    Lout_feat(i,8) = {mean(gabor_extract)}; %average gabor
    Lout_feat(i,9) = {[mean(X_coord(L_out_mask)), mean(Y_coord(L_out_mask))]}; %average distance
end

X_all = get_feature(Lin_feat, Lout_feat, cov_feat,num_samp_Lin, u_Lout, offset);



if ~(isempty(chi_match) || isempty(chi_distance)) % Training
    %pos = 0.075;
    %neg = 0.075;
    %{
    max_count = -1;
    for i = 1:num_samp_Lin
        if(max_count<sum(clustered_L_in(:)==i)), max_count=sum(clustered_L_in(:)==i); end
    end
    %}
    count = 1;
    for i = 1:num_samp_Lin
        pos = max(1.0, chi_distance(i,1)+0.5); % max(1.5, chi_distance(i,1)+0.5); %max(0.05, chi_distance(i,1)+0.01);
        neg = min(pos+1.5, chi_distance(i,end)-0.5); %min(pos+1.5, chi_distance(i,end)-0.5); %min(pos+0.2, chi_distance(i,end)-0.01);
        %if(sum(clustered_L_in(:)==i)<min(3e2, max_count)), continue; end
        p_idx = chi_match(i, chi_distance(i,:)<pos);    
        n_idx = chi_match(i, chi_distance(i,:)>neg);
        pn_idx = [p_idx, n_idx]; pn = [ones(length(p_idx),1); zeros(length(n_idx), 1)];
        for j = 1:size(pn_idx,2)
            %X(i, pn_idx(j), 1:size(X_all,3)) = X_all(i, pn_idx(j), :);
            %Y(i, pn_idx(j), 1) = pn(j);            
            X(count, 1:size(X_all,3)) = X_all(i, pn_idx(j), :);
            Y(count, 1) = pn(j);
            count = count+1;
        end
    end
    %Y = reshape(Y,num_samp_Lin* max(u_Lout(:)),1);
    %X = reshape(X, num_samp_Lin* max(u_Lout(:)), []);
else
    %Testing
    %X = zeros(num_samp_Lin,max(u_Lout(:)), size(X_all,3));
    Y = zeros(max(u_Lout(:))*num_samp_Lin, 1);
    X = reshape(X_all, num_samp_Lin, max(u_Lout(:)), []);
end




end

function [angles, magnit] = gradient_profile(gray_img)

% Create the operators for computing image derivative at every pixel.
hx = 0.5*[-1,0,1];
hy = hx';

% Compute the derivative in the x and y direction for every pixel.
dx = filter2(hx, padarray(double(gray_img), [1,1]));
dy = filter2(hy, padarray(double(gray_img), [1,1]));

% Remove the 1 pixel border.
dx = dx(2 : (size(dx, 1) - 1), 2 : (size(dx, 2) - 1));
dy = dy(2 : (size(dy, 1) - 1), 2 : (size(dy, 2) - 1)); 

% Convert the gradient vectors to polar coordinates (angle and magnitude).
angles = atan2(dy, dx);
magnit = ((dy.^2) + (dx.^2)).^.5;
end


function [histw] = histcounts_w(vv, ww, edges)
  histw = zeros(1,length(edges));
  for i=1:length(vv)
    ind = find(edges < vv(i), 1, 'last' );
    if ~isempty(ind)
      histw(ind) = histw(ind) + ww(i);
    end
  end
  histw = histw./sum(histw);
end
