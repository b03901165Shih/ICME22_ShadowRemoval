% Use gt image to find positive and negative matches
function [ chi_match, chi_distance] = get_region_matching( gt, clustered_L_in, clustered_L_out )

%Construct color histogram 

u_Lin = unique(clustered_L_in);
u_Lout = unique(clustered_L_out);

I_lab = rgb2lab(gt);
I_lab(:,:,1) = I_lab(:,:,1)/100;
for c = 2:3, I_lab(:,:,c) = (I_lab(:,:,c)+100)/(200); end

r = gt(:,:,1); g = gt(:,:,2); b = gt(:,:,3);
l_s = I_lab(:,:,1); a_s = I_lab(:,:,2); b_s = I_lab(:,:,3);

edges = 0:0.02:1;
%nbins=50;
Hist_Lin = cell(max(u_Lin),2);
Hist_Lout = cell(max(u_Lout),2);
    
for i = 1:max(u_Lin(:))
    select_map = clustered_L_in==i;
    Hist_Lin(i, 1) = {cat(1, ...
        histcounts(l_s(select_map), edges,'Normalization', 'probability'),  ...
        histcounts(a_s(select_map), edges, 'Normalization', 'probability'),  ...
        histcounts(b_s(select_map), edges, 'Normalization', 'probability'))+eps};
    %Hist_Lin(i, 2) = {mean(reshape(DAISY_f(repmat(select_map, [1 1 g_dim])), [], g_dim))};
end
for i = 1:max(u_Lout(:))
    select_map = clustered_L_out==i;
    Hist_Lout(i, 1) = {cat(1, ...
        histcounts(l_s(select_map), edges,'Normalization', 'probability'),  ...
        histcounts(a_s(select_map), edges, 'Normalization', 'probability'),  ...
        histcounts(b_s(select_map), edges, 'Normalization', 'probability'))+eps};
	%Hist_Lout(i, 2) = {mean(reshape(DAISY_f(repmat(select_map, [1 1 g_dim])), [], g_dim))};
end

chi_distance = zeros(max(u_Lin(:)), max(u_Lout(:)));
%gabor_distance = zeros(max(u_Lin(:)), max(u_Lout(:)));
for i = 1:max(u_Lin(:))
    for j = 1:max(u_Lout(:))
        D=0;
        for c = 1:size(Hist_Lin{i, 1},1)
            D = D+pdist2_mod( Hist_Lin{i, 1}(c,:), Hist_Lout{j, 1}(c,:),  'emd');
        end
        chi_distance(i,j) = D/size(Hist_Lin{i, 1},1); %(D/size(Hist_Lin{i, 1},1)+sqrt(sum( (Hist_Lin{i,2}- Hist_Lout{j,2}).^2)))/2;
        %gabor_distance(i,j) = sqrt(sum( (Hist_Lin{i,2}- Hist_Lout{j,2}).^2)); %
    end
end

chi_match = zeros(max(u_Lin(:)), max(u_Lout(:)));
for i = 1:max(u_Lin(:))
    [sort_val, sort_i]=sort(chi_distance(i,:));
    chi_match(i,:) = sort_i;
    chi_distance(i,:) = sort_val;
end

end

