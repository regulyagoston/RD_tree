

function sq_l = get_np_split_criterion(Y, X, bw_j, c, z_split_id, fuzzy, cl, optTree, N_tr, N_est, z_split_id_est )

if all( isnan( fuzzy ) )
    fuzzy = [];
else
    fuzzy = fuzzy( z_split_id );
end
if all( isnan( cl ) )
    cl = [];
else
    cl = cl( z_split_id );
end

out_l = rdrobust( Y( z_split_id ) , X( z_split_id ) , c , bw_j, fuzzy , cl ,...
                optTree.orderPolinomial, optTree.varEstimator, optTree );
n_j_tr_l = sum( z_split_id );
p_j_tr_l = n_j_tr_l/N_tr;
n_j_est_l = sum( z_split_id_est );
p_j_est_l = n_j_est_l/N_est;
% Expected value of SQUARED treatments
e_tau_sq_l = -( out_l.tau_bc.^2 - out_l.tau_V_rb * p_j_tr_l );
% Bias
bias_sq_l = bw_j.h_l.^2 .* (out_l.bias_r - out_l.bias_l).^2 .* p_j_est_l;
% Variance
var_j_l = out_l.tau_V_cl*n_j_tr_l/N_est;%*optTree.bw/optTree.bw;
sq_l = e_tau_sq_l + bias_sq_l + var_j_l;





end