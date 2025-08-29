%% Predict the leaves of a tree

function obj_tree_pred = predictTree( obj_tree_pred , obj_tree_tr , obj_sample , obj_optCART )

lfID = findleaves( obj_tree_pred );
nL   = numel( lfID );

for i = 1 : nL
    obj_tree_pred.Node{ lfID( i ) } = predictNode( obj_tree_pred.Node{ lfID( i ) } ,...
                                                   obj_tree_tr.Node{ lfID( i ) } , ...
                                                   obj_sample , obj_optCART );
end


end