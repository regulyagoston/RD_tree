

function [ est_tree , obj_tree , opt_beta , s_oOS, opt_bw, opt_bw_b ] = runTree_uni( obj_sample , optTree, cv_type, num_bw )


% Cross-validation to find optimal complexity parameter (beta) and
% bandwidth
optTree = copyWithChange(optTree,'num_bw', num_bw );
if strcmp( cv_type, 'fmincon' )
    [ opt_beta , opt_bw, s_oOS, large_tree, rho ] = cross_validate_bws( obj_sample , optTree );
    optTree.bw = opt_bw;
    optTree.bw_b = opt_bw .* rho;
    optTree.num_bw = num_bw;
elseif strcmp( cv_type, 'grid' )
    [ opt_beta , opt_bw, opt_bw_b, s_oOS ] = cross_validate_bws_grid( obj_sample , optTree, [] );
    optTree = copyWithChange(optTree,'bw', opt_bw);
    optTree = copyWithChange(optTree,'bw_b', opt_bw_b);
    large_tree = growRDDtree( obj_sample , optTree );
else
    error('Invalid input for cv_type!')
end
obj_tree = pruning( large_tree , obj_sample , optTree , opt_beta );
%% Estimate proper CATE
est_tree = estimate_tree_est( obj_tree , obj_sample , optTree );

end