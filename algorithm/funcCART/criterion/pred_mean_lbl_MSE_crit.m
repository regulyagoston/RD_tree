%% Prediction criterion using MSE

function sQ = pred_mean_lbl_MSE_crit( leaf_nodes , obj_sample , obj_optCART )

out_struct = collectInputs( leaf_nodes , true , 'Y_hat_i' );

Y_hat_i = out_struct;

sQ = mean( ( obj_sample.Y_tr - Y_hat_i ) .^ 2 );

end