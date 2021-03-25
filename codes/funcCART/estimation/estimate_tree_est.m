%% Estimate the parameters of interests for the tree using the estimation sample



function est_tree = estimate_tree_est( old_tree , obj_sample , obj_optCART )

cleared_tree = clearTree( old_tree );
est_tree = updateTree( cleared_tree , old_tree , obj_sample , obj_optCART , false );

end