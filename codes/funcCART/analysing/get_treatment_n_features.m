%% Get treatment effect along one feature


function [ tau , tau_se , Z ] = get_treatment_n_features( tree , sample )


lID = findleaves( tree );
nL = numel( lID );

tau = sample.n_est;
tau_se = sample.n_est;
Z = NaN( size( sample.Z_est ) );

for i = 1 : nL
    leaf_log = tree.Node{ lID( i ) }.logID_est;
    tau( leaf_log ) = tree.Node{ lID( i ) }.est_struct.tau_i;
    tau_se( leaf_log ) = tree.Node{ lID( i ) }.est_struct.tau_se_j;
    Z( leaf_log , : ) = sample.Z_est( leaf_log , : );
end

end