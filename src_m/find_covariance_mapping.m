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

% RGB version
%%{
function [ A, B] = find_covariance_mapping( ref_vec, target_vec)

% numbers of data
N1 = size(ref_vec, 1);
N2 = size(target_vec, 1);

%Matching mean only
mean_ref = mean(ref_vec, 1);
mean_target = mean(target_vec, 1);

std_ref = std(ref_vec, 1);
std_target = std(target_vec, 1);

A = diag(std_ref./std_target);
B = mean_ref-mean_target.*(std_ref./std_target);
end
%}

