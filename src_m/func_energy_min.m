function [ f_t , foreground_remapped, E] = func_energy_min(params, foreground, remapped_f, clustered_L_in, clustered_L_out, Y )

match_idx = zeros(size(Y));
match_score = zeros(size(Y));
for i = 1:size(Y,1)
    [match_score(i,:), match_idx(i,:)]=sort(Y(i,:), 'descend');
end

%Global mapping function
%{
RGB_shadow_all = reshape(foreground(repmat(clustered_L_in~=0,[1 1 3])),[],3);
RGB_lit_all = reshape(foreground(repmat(ismember(clustered_L_out, unique(match_idx(match_score>=0.5))),[1 1 3])),[],3);
[A_global, B_global] = find_covariance_mapping( RGB_lit_all, RGB_shadow_all); 
%}

RGB_shadow = cell(size(Y, 1), 1);
RGB_lit = cell(size(Y, 2), 1);
for i = 1:size(Y, 1)
    RGB_shadow(i) = {reshape(foreground(repmat(clustered_L_in==i,[1 1 3])),[],3)};
end
for i = 1:size(Y, 2)
    RGB_lit(i) = {reshape(foreground(repmat(clustered_L_out==i,[1 1 3])),[],3)};
end


% Similar weight
phi = zeros(size(Y,1));
[ chi_match, chi_distance ] = get_region_matching( foreground, clustered_L_in, clustered_L_in );
for i = 1:size(Y, 1)
    [~, sort_i]=sort(chi_match(i,:));
    sorted_dis = chi_distance(i,:);
    phi(i,:) = exp(-sorted_dis(sort_i)/1.0); %/0.5
end

% 2-dim: A/ B
f_init = cell(size(Y, 1), 2);
f_p = cell(size(Y, 1), 2);
for i  = 1:size(Y, 1)
    RGB_lit_part = reshape(remapped_f(repmat(clustered_L_in==i,[1 1 3])),[],3);
    [A, B] = find_covariance_mapping( RGB_lit_part, RGB_shadow{i}); %RGB_lit{Klist(k)}
    f_init(i,1) = {A};
    f_init(i,2) = {B};
    if(match_score(i,1)<0.5)
        continue;
    end
    %
    Klist = match_idx(i,match_score(i,:)>=0.5);%match_idx(i,v{i}~=0)
    A_k = cell(length(Klist), 1); B_k = cell(length(Klist), 1);
    for k = 1:length(Klist)
        [A, B] = find_covariance_mapping( RGB_lit{Klist(k)}, RGB_shadow{i}); %RGB_lit{Klist(k)}
        A_k(k) = {A};
        B_k(k) = {B};
    end
    f_p(i,1) = {A_k};
    f_p(i,2) = {B_k};
end

% RGB_shadow: cell(size(Y, 1), 1);
% f_init: cell(size(Y, 1), 2)
% f_p: cell(size(Y, 1), 2) of cell(length(Klist), 1)
% phi: zeros(size(Y,1), size(Y,1));
%%  optimization

alpha_p = 1;
gamma_p = 0.0005;
beta_p = gamma_p;
beta_max = gamma_p*2^10;
epsilon_p = 1e-2;
options = optimset('Display','none');

w = cell(size(Y,1), 1);
f_t = f_init;
% Initialize u, w, and u_similar
for i = 1:size(Y,1)
    w(i) = {1/length(f_p{i})*ones(length(f_p{i}),1)};
end

% Average # of pixels for a cluster
norm_factor = 1/(sum(clustered_L_in(:)~=0)/max(unique(clustered_L_in)));

sum_energy_l = [];
match_energy_l = [];
sim_energy_l = [];
l0_energy_l = [];

sum_energy_all_old = inf;
sum_energy_old = inf;
diff_ratio_all = 1;
max_it = 50;
count_it_out = -1;
while(diff_ratio_all > epsilon_p && count_it_out < max_it)
    count_it_out = count_it_out+1;
    % f step
    diff_ratio = 1;
    f_t_old = f_t;
    count_it_in = -1;
    while(diff_ratio>epsilon_p && count_it_in < max_it)
        count_it_in = count_it_in+1;
        % A update
        for i = 1:size(Y, 1)
            lhs_mat = zeros(3,3); rhs_mat = zeros(3,3);
            U_p = RGB_shadow{i}; B_p = f_t{i,2};
            one_p = ones(size(U_p,1),1);
            for k = 1:length(f_p{i}) 
                A_pk = f_p{i,1}{k}; B_pk = f_p{i,2}{k}; 
                lhs_mat = lhs_mat+norm_factor*w{i}(k)*(U_p'*U_p);%+eps*eye(3);    
                rhs_mat = rhs_mat+norm_factor*w{i}(k)*U_p'*bsxfun(@minus, U_p*A_pk, (B_p-B_pk));
            end
            for j = 1:size(Y, 1)
                U_p = RGB_shadow{j}; B_p = f_t{j,2};
                one_p = ones(size(U_p,1),1);
                A_psim = f_t_old{j,1}; B_psim = f_t_old{j,2}; 
                lhs_mat = lhs_mat+norm_factor*alpha_p*phi(i,j)*(U_p'*U_p);%+eps*eye(3);    
                rhs_mat = rhs_mat+norm_factor*alpha_p*phi(i,j)*U_p'*bsxfun(@minus, U_p*A_psim, (B_p-B_psim));
            end
            f_t(i,1) = {lhs_mat\rhs_mat};
        end
        % B update
        for i = 1:size(Y, 1)
            lhs_mat = 0; rhs_mat = zeros(1,3);
            U_p = RGB_shadow{i}; A_p = f_t{i,1};
            one_p = ones(size(U_p,1),1);
            for k = 1:length(f_p{i}) 
                A_pk = f_p{i,1}{k}; B_pk = f_p{i,2}{k}; 
                lhs_mat = lhs_mat+norm_factor*w{i}(k)*size(U_p,1);%(one_p'*one_p);
                rhs_mat = rhs_mat+norm_factor*w{i}(k)*sum(bsxfun(@plus, U_p*(A_pk-A_p), (B_pk)));
            end
            for j = 1:size(Y, 1)
                U_p = RGB_shadow{j}; A_p = f_t{j,1};
                one_p=ones(size(U_p,1),1);
                A_psim = f_t_old{j,1}; B_psim = f_t_old{j,2}; 
                lhs_mat = lhs_mat+norm_factor*alpha_p*phi(i,j)*size(U_p,1);%(one_p'*one_p);
                rhs_mat = rhs_mat+norm_factor*alpha_p*phi(i,j)*sum(bsxfun(@plus, U_p*(A_psim-A_p), (B_psim)));
            end
            f_t(i,2) = {rhs_mat/lhs_mat};
        end

        sum_energy = 0;
        for i = 1:size(Y, 1)
            U_p = RGB_shadow{i}; 
            one_p = ones(size(U_p,1),1);
            for k = 1:length(f_p{i}) 
                sum_energy = sum_energy+norm_factor*1/2*w{i}(k)*norm(U_p*(f_t{i,1}-f_p{i,1}{k})+one_p*(f_t{i,2}-f_p{i,2}{k}),'fro')^2;
            end
            for j = 1:size(Y, 1)
                U_p = RGB_shadow{j}; 
                one_p = ones(size(U_p,1),1);
                sum_energy = sum_energy+norm_factor*1/2*alpha_p*phi(i,j)*norm(U_p*(f_t{i,1}-f_t_old{j,1})+one_p*(f_t{i,2}-f_t_old{j,2}),'fro')^2;
            end
        end
        diff_ratio = abs(sum_energy_old-sum_energy)/sum_energy;
        sum_energy_old = sum_energy;
        if(params.verbose>1)
            fprintf('\tEnergy : %2.5f | Energy diff :%2.5f (stopping at %2.5f)\n', sum_energy, diff_ratio, epsilon_p);
        end
    end

    %% w-subproblem
    % initilize v
    v = w;
    for i = 1:size(Y,1)
        if(isempty(w)), continue; end
        len_w = length(w{i});
        beta_curr = beta_p*len_w^2;
        while(beta_curr<(beta_max*len_w^2))
            % w sub-sub-problem
            g = zeros(len_w, 1);
            U_p = RGB_shadow{i};
            lit_pixel_i = bsxfun(@plus, U_p*f_t{i,1}, f_t{i,2}); 
            for k = 1:len_w
                g(k) = norm_factor*norm(lit_pixel_i-bsxfun(@plus, RGB_shadow{i}*f_p{i,1}{k}, f_p{i,2}{k}),'fro')^2;
            end
            f_quad = 1/2*g-beta_curr*v{i};
            H_quad = beta_curr*eye(len_w); 
            Aeq = ones(1,len_w); beq = 1; lb=zeros(len_w,1); 
            w(i) = {quadprog(H_quad,f_quad,[],[],Aeq,beq,lb,[],w{i},options)};
            % v sub-sub-problem
            %mask_v = (w{i}.^2)>(2*gamma_p/beta_curr);
            %v(i) = {w{i}.*mask_v};
            v(i) = {wthresh(w{i},'h',sqrt(2*gamma_p/beta_curr))};
            % increase beta
            beta_curr = beta_curr*2;
        end
    end
    
    match_energy = 0; sim_energy = 0; l0_energy = 0;
    for i = 1:size(Y, 1)
        U_p = RGB_shadow{i}; 
        one_p = ones(size(U_p,1),1);
        for k = 1:length(f_p{i}) 
            match_energy = match_energy+norm_factor*1/2*w{i}(k)*norm(U_p*(f_t{i,1}-f_p{i,1}{k})+one_p*(f_t{i,2}-f_p{i,2}{k}),'fro')^2;
        end
        for j = 1:size(Y, 1)
            U_p = RGB_shadow{j}; 
            one_p = ones(size(U_p,1),1);
            sim_energy = sim_energy+norm_factor*1/2*alpha_p*phi(i,j)*norm(U_p*(f_t{i,1}-f_t_old{j,1})+one_p*(f_t{i,2}-f_t_old{j,2}),'fro')^2;
        end
        if(~isempty(v{i})), l0_energy = l0_energy+gamma_p*sum(v{i}~=0); end
    end
    sum_energy_all = match_energy+sim_energy+l0_energy;
    
    sum_energy_l = [sum_energy_l, sum_energy_all];
    match_energy_l = [match_energy_l, match_energy];
    sim_energy_l = [sim_energy_l, sim_energy];
    l0_energy_l = [l0_energy_l, l0_energy];

    %%
    diff_ratio_all = abs(sum_energy_all_old-sum_energy_all)/sum_energy_all;
    sum_energy_all_old = sum_energy_all;
    if(params.verbose>1)
        fprintf('All Energy : %2.5f | Energy diff :%2.5f (stopping at %2.5f) | Match energy: %2.5f | Sim energy: %2.5f | L0 energy: %2.5f\n', ...
        sum_energy_all, diff_ratio_all, epsilon_p, match_energy,sim_energy,l0_energy);
    end
end

E = sum_energy_all;

%{
figure, plot(sum_energy_l), xlabel('Iteration'), ylabel('Energy'), grid on, hold on;
plot(match_energy_l), hold on;
plot(sim_energy_l), hold on;
plot(l0_energy_l), hold on;
legend('Total', 'Match','Similarity', 'L0');
%}

% Reconstruct f_t from only the matching
%{
f_t_old = f_t;
for i  = 1:size(Y, 1)
    if(match_score(i,1)<0.5)
        continue;
    end
    Klist = match_idx(i, v{i}~=0);
	outer_mask = imclose(ismember(clustered_L_out,Klist), ones(10));%match_score(i,:)>=0.5
	RGB_lit_part = reshape(foreground(repmat(outer_mask, [1,1,3])), [], 3);
    [A, B] = find_covariance_mapping( RGB_lit_part, RGB_shadow{i}); %RGB_lit{Klist(k)}
    f_t(i,1) = {A};
    f_t(i,2) = {B};
end
%}

%% Visualize f_t
foreground_remapped = foreground;
[yidx, xidx]=meshgrid(0:size(foreground, 2)-1, 0:size(foreground, 1)-1);
for i = 1:size(Y, 1)
    L_in_mask = (clustered_L_in==i);
    remap_c = bsxfun(@plus, RGB_shadow{i}*(f_t{i,1}), f_t{i,2});
    for ch = 1:3
        rgb = foreground_remapped(:,:,ch);
        rgb(xidx(L_in_mask)+yidx(L_in_mask)*size(foreground, 1)+1) = remap_c(:,ch);
        foreground_remapped(:,:,ch) = rgb;
    end
end

%{

for i = 1:min(10, size(Y, 1))
    idx = match_idx(i,v{i}~=0);
    figure, imshow([foreground.*repmat(clustered_L_in==i,[1,1,3])+(0.2+0.2*foreground).*repmat(clustered_L_in~=i,[1,1,3]), ...
    foreground.*repmat(ismember(clustered_L_out, idx),[1,1,3])+(0.2+0.2*foreground).*repmat(~ismember(clustered_L_out, idx),[1,1,3])])
    %imwrite([foreground.*repmat(clustered_L_in==i,[1,1,3])+(0.2+0.2*foreground).*repmat(clustered_L_in~=i,[1,1,3]), ...
    %foreground.*repmat(ismember(clustered_L_out, idx),[1,1,3])+(0.2+0.2*foreground).*repmat(~ismember(clustered_L_out, idx),[1,1,3])], ['1152-after', int2str(i),'.png']);
end
%}
%% Finding the mapping function
%[ A, B] = find_covariance_mapping( ref_vec, target_vec)
%{

mapf_init = cell(size(Y, 1), 1);
for i  = 1:size(Y, 1)
    Klist = match_idx(i,match_score(i,:)>=0.5);
    weight_A = zeros(3,3); weight_B = zeros(1,3);
	%outer_mask = imclose(ismember(clustered_L_out,match_idx(i,match_score(i,:)>=0.5)), ones(10));%match_score(i,:)>=0.5
	%RGB_lit_part = reshape(foreground(repmat(outer_mask, [1,1,3])), [], 3);
    for k = 1:length(Klist)
        [A, B] = find_covariance_mapping( RGB_lit_part, RGB_shadow{i}); %RGB_lit{Klist(k)}
        weight_A = weight_A+A/length(Klist);
        weight_B = weight_B+B/length(Klist);
    end
    func.A = weight_A;
    func.B = weight_B;
    mapf_init(i) = {func};
end


%[mapf] = function_energy_minimization(RGB_shadow, similar_pk, mapf_init);
mapf = mapf_init;

foreground_remapped = foreground;
[yidx, xidx]=meshgrid(0:size(foreground, 2)-1, 0:size(foreground, 1)-1);
for i = 1:size(Y, 1)
    L_in_mask = (clustered_L_in==i);
    remap_c = bsxfun(@plus, RGB_shadow{i}*(mapf{i}.A), mapf{i}.B);
    for ch = 1:3
        rgb = foreground_remapped(:,:,ch);
        rgb(xidx(L_in_mask)+yidx(L_in_mask)*size(foreground, 1)+1) = remap_c(:,ch);
        foreground_remapped(:,:,ch) = rgb;
    end
end
figure, imshow(foreground_remapped)
%}

end