%% Prediction criterion using MSE only for splitting candidates

function sQ = pred_mean_lbl_MSE_SPLIT_crit( split_nodes , obj_sample , obj_optCART )

% Take all the training observations
Y_tr = obj_sample.Y_tr;
% Create vector for estimated values
Y_hat_i = NaN( size( Y_tr ) );
% Get the id-s for the splitting nodes
id_1 = split_nodes(1).logID_tr;
id_2 = split_nodes(2).logID_tr;
% Id for the both
id   = id_1 | id_2;
% Put the estimated values into the vector
Y_hat_i( id_1 ) = split_nodes( 1 ).est_struct.Y_hat_i;
Y_hat_i( id_2 ) = split_nodes( 2 ).est_struct.Y_hat_i;

% Calculate MSE only on the candidates
sQ = mean( ( Y_tr( id ) - Y_hat_i( id ) ) .^ 2 );

end