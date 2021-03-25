%% Grow, cross-validate, prune and estimate RDD tree


function [ est_tree , obj_tree , opt_beta , s_oOS , beta , oOS ] = runTree( obj_sample , optTree )


[ obj_tree , stopWhy ] = growRDDtree( obj_sample , optTree );
% Cross-validation
[ opt_beta , s_oOS , beta , oOS ] = cross_validate( obj_tree , obj_sample , optTree );
% Prune with optimal complexity parameter
finalTree_tr = pruning( obj_tree , obj_sample , optTree , opt_beta );
% Do the estimation on the estimation sample
est_tree = estimate_tree_est( finalTree_tr , obj_sample , optTree );

end