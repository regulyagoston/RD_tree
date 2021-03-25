%% Clear the tree's nodes

function obj_tree = clearTree( obj_tree )

nN     = nnodes( obj_tree );
% Make sure it is not an empty tree
if ~isempty( obj_tree.Node{ 1 } )
    for i = 1 : nN
        % Clear all the estimated parameters for the node
        obj_tree.Node{ i } = clearNode( obj_tree.Node{ i } );
    end
else
    error('tree:clearTree:incorrectType',...
          'Can not clear the tree, while it is empty!')
end

end