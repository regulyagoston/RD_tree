%% Update the leaves of a tree

function obj_tree = updateTree( obj_tree , old_tree , obj_sample , obj_optCART , tr )

nL   = nnodes( obj_tree );

for i = 1 : nL
    try
        obj_tree.Node{ i } = updateNode( obj_tree.Node{ i } , obj_sample , obj_optCART , tr );
    catch
        % warning('There is no enough observation in the leaf, predicting values instead of estimating!');
        obj_tree.Node{ i } = predictNode( obj_tree.Node{ i } , old_tree.Node{ i } , obj_sample , obj_optCART , tr );
    end
end


end