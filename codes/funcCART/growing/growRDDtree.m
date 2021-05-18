%% Growing large tree
% Agoston Reguly (2020)
%
% Input:
%   data    - samples object
%   optTree - optCART object
%
% Output:
%   obj - tree object
%   stopReason - string
%       why the algorithm stopped

function [ obj , stopReason ] = growRDDtree( data , optTree )

%% Creating the Tree
obj = tree;
% Create the root node, where all observations are present
rootNode = nodeProp( [] , [] , [] , 0 , true( 1 , data.K ) );
% Update the node -> according the the criterion. Growing trees are always
%    done on the training sample
rootNode = updateNode( rootNode , data , optTree , true );
% Calculate the ios for the rootnode
% ios_0 = optTree.crit_ios( rootNode , data , optTree );
% Add Properties and calculate the criterion
obj = addnode( obj , 0 , rootNode );
% Set parent to 1 to create children
parent = 1;

%% Start to grow the tree
% Notes:
% Only split to two parts at each iteration
% Once a split is made it is fixed. In next iteration only searches along
%   the already splitted sample

% Stop if any is satisfied
if any( [ optTree.maxLevel == 1 , optTree.maxNodes == 1 , optTree.maxLeaves == 1 , ...
        ( 2 * optTree.minObs ) > obj.Node{ 1 }.n_j_tr ] )
    stopReason = 'No more valid nodes';
    return;
end

numNode_j = 1;
numLeaves = 1;
IDnode    = 1;
IDsib     = 1;
for j = 1 : optTree.maxIterGrow
    
    %% Get potential variables to split
    % 1st - select those features which are a valid split on that node
    [ id_valid , n_id ] = get_valid_features( obj.Node{ parent } );
    % 2nd - create n_id number of nodes for the potential features
    cand_bS_f = repmat( nodeProp , [ n_id , 2 ] );
    % 3rd - assume that all potential feature is valid
    cand_valid = true( 1 , n_id );
    cand_ios   = NaN( 1 , n_id );
    % 4th - Find the best split for each candidate feature
    try
        if ~optTree.paralell
            error('This is only a shortcut if no paralell processing is used!')
        end
        parfor k = 1 : n_id
            [ cand_bS_f( k , : ) , cand_valid( k ) , cand_ios( k ) ] = ...
                    best_split_by_feature( obj , parent , id_valid( k ) , data , optTree );
        end
    catch
       for k = 1 : n_id
            [ cand_bS_f( k , : ) , cand_valid( k ) , cand_ios( k ) ] = ...
                    best_split_by_feature( obj , parent , id_valid( k ) , data , optTree );
        end 
    end
        
    % 5th - Save the result to the parent's node:
    % If failed, why...
    %obj.Node{ parent }.feature_chck( id_valid ) = cand_valid;
    % Mark that this node has been investigated
    obj.Node{ parent }.Checked = true;
    
    %% Select the best feature among the candidates
    if any( cand_valid )
        % Select the node with best split criterion
        bestCrit = min( cand_ios( cand_valid ) );
        bestNode = cand_bS_f( cand_ios == bestCrit , : );
        % Add Nodes to the tree
        [ obj , IDnode ] = addnode( obj , parent , bestNode( 1 , 1 ) );
        [ obj , IDsib  ] = addnode( obj , parent , bestNode( 1 , 2 ) );
    end
    
    %% Check for stopping the algorithm
    % Check for depth of the tree. If reached the maximum level set it
    % a terminal node
    if ~isinf( optTree.maxLevel )
        dt = depthtree( obj );
        node_level = dt.Node{ IDnode };
        if node_level == optTree.maxLevel
            obj.Node{ IDnode }.Checked = true;
            obj.Node{ IDnode }.Failed  = 3;
            obj.Node{ IDsib }.Checked = true;
            obj.Node{ IDsib }.Failed  = 3;
        end
    end
    % Too much nodes -> set it as a terminal node 
    if optTree.maxNodes <= numNode_j
        obj.Node{ IDnode }.Checked = true;
        obj.Node{ IDnode }.Failed  = 4;
        stopReason = 'Tree-size reached the maximal number of nodes, tree growing stopped!';
        %disp( 'Tree-size reached the maximal number of nodes, tree growing stopped!' )
        return;
    end
    % Too much leaves -> set it as a terminal node
    if optTree.maxLeaves <= numLeaves
        obj.Node{ IDnode }.Checked = true;
        obj.Node{ IDnode }.Failed  = 5;
        stopReason = 'Tree-size reached the maximal number of leaves, tree growing stopped!';
        %disp( 'Tree-size reached the maximal number of leaves, tree growing stopped!' )
        return;
    end
    
    % Not enough observations in the leaf -> set it as a terminal node
    if obj.Node{ IDnode }.n_j_tr < 2 * optTree.minObs
        obj.Node{ IDnode }.Checked = true;
        obj.Node{ IDnode }.Failed  = 6;
    end
    if obj.Node{ IDsib }.n_j_tr < 2 * optTree.minObs
        obj.Node{ IDsib }.Checked = true;
        obj.Node{ IDsib }.Failed  = 6;
    end
    
    
    % Set the next candidate parent:
    %   next candidate split is the next valid child
    %   this means the algorithm searches through levels first and than
    %   creates a new 
    parent = parent + 1;
    
    % Update the number of nodes and leaves
    numNode_j = nnodes( obj );
    numLeaves = numel( findleaves( obj ) );
    
    %% Check and get valid nodes
    % No more valid nodes
    if parent > numNode_j
        stopReason = 'No more valid nodes';
        break;
    end
    
    % Get valid parent node
    while obj.Node{ parent }.Checked
        parent = parent + 1;
        if parent > numNode_j
            % No more valid nodes
            stopReason = 'No more valid nodes';
            return;
        end
    end
    
end

end