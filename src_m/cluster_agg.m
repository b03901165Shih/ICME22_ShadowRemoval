function [ clustered_L_in, clustered_L_out ] = cluster_agg( params, foreground, L_in, L_out)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    clustered_L_in = L_in;
    clustered_L_out = L_out;
    if(params.verbose>0)
        fprintf('Clusters Aggregation\n');
    end
    if(params.aggregate_cluster)
        threshold = 0.10; 
        while(max(unique(clustered_L_in))>30)
            [ clustered_L_in ] = clustering_using_aggtree( L_in, foreground, threshold);
            threshold = threshold+0.025;
        end    
        threshold = 0.10; 
        while(max(unique(clustered_L_out))>30)
            [ clustered_L_out ] = clustering_using_aggtree( L_out, foreground, threshold);
            threshold = threshold+0.025;
        end
    end    
    if(params.verbose>0)
        fprintf('In/out number of clusters: %d, %d\n', max(unique(clustered_L_in)), max(unique(clustered_L_out)));    
    end

end

function [ clustered_L ] = clustering_using_aggtree( L, foreground, threshold)

    u_L = unique(L);


    % center pixel coord
    [Y, X] = meshgrid(1/size(foreground, 2):1/size(foreground, 2):1, 1/size(foreground,1):1/size(foreground,1):1);
    center_sp = zeros(max(u_L(:)), 2);
    for i = 1:max(u_L(:))
        sp_mask = L ==i;
        center_sp(i,:) = [mean(X(sp_mask)), mean(Y(sp_mask))];
    end

    %%
    edges = 0:0.02:1;
    dim = length(edges)-1;
    feat_L = zeros(max(u_L(:)), 3*dim);

    %I_lab = rgb2ycbcr(foreground);
    I_lab = rgb2lab(foreground);
    I_lab(:,:,1) = I_lab(:,:,1)/100;
    for c = 2:3, I_lab(:,:,c) = (I_lab(:,:,c)+100)/(200); end

    %dist_ch = @(ZI, ZJ, ch) (1/max_ch)*pdist2_mod(ZI(:,dim*(ch-1)+1:dim*ch), ZJ(:,dim*(ch-1)+1:dim*ch), 'emd')';
    %distfun = @(ZI,ZJ) (dist_ch(ZI,ZJ,1)+dist_ch(ZI,ZJ,2)+dist_ch(ZI,ZJ,3))/3;

    for i = 1:max(u_L(:))
        pixels = reshape(I_lab(repmat(L==i, [1,1,3])),[],3);
        for ch = 1:3
            feat_L(i, dim*(ch-1)+1:dim*ch) = histcounts(pixels(:,ch), edges,'Normalization', 'probability');
        end
    end
    feat_L = cat(2, feat_L, center_sp);


    Z = linkage(feat_L,'weighted', {@distfun}); %complete 'single' 'weighted'
    C = cluster(Z,'cutoff',threshold,'Criterion','distance');%'MaxClust',10);

    %% update L
    clustered_L = L;
    for n = 1:max(L(:))
        clustered_L(L==n) = C(n);
    end

end

function D = distfun(ZI,ZJ)
    D = zeros(size(ZI,1), size(ZJ,1));
    max_ch =3; dim = 50;
    for ch = 1:max_ch
        D = D+(1/max_ch)*pdist2_mod(ZI(:,dim*(ch-1)+1:dim*ch), ZJ(:,dim*(ch-1)+1:dim*ch), 'emd');
    end
    spatial = pdist2(ZI(:,end-1:end), ZJ(:,end-1:end)); %1-exp(-pdist2(ZI(:,end-1:end), ZJ(:,end-1:end)));
    D = (D.*spatial)';
end

