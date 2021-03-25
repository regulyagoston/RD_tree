%% Update the leaves of a tree

function obj_tree = updateTree_leaves( obj_tree , obj_sample , obj_optCART , tr )

lfID = findleaves( obj_tree );
nL   = numel( lfID );

for i = 1 : nL
    obj_tree.Node{ lfID( i ) } = updateNode( obj_tree.Node{ lfID( i ) } , obj_sample , obj_optCART , tr );
end


end