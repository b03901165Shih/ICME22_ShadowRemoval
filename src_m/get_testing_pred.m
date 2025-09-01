function [ Y ] = get_testing_pred( foreground, clustered_L_in, clustered_L_out, Mdl)

u_Lin = unique(clustered_L_in);
u_Lout = unique(clustered_L_out);

[ X, ~ ] = get_training_feature( foreground, [], [], clustered_L_in, clustered_L_out );

reshape_X = reshape(X, [], size(X, 3));
reshape_Y = predict(Mdl, reshape_X);
Y = reshape(reshape_Y, max(u_Lin(:)), max(u_Lout(:)));

end

