%% Prune the tree
% If alpha is given, prune back with alpha
%   if alpha is not given, find the pruning values (alphas)
%
% Input:
%   obj_tree - tree object
%       contains the tree to be pruned
%   obj_sample - sample object
%       contains the data
%   obj_optCART - optCART object
%       contains options
%   alpha0 - optional -> scalar
%       If alpha is given, prune back with alpha
%       if alpha is not given, find the pruning values (alphas)
%
% Output:
%   prTrees - vector/scalar tree object
%       pruned tree(s)
%   iOS - vector/scalar
%       gives the in-sample criterion's value - relative values to the
%       fully grown tree, which has value of 0
%   alpha - vector/scalar
%       gives the optimal alpha values


function [ prTrees , iOS , alpha ] = pruning( obj_tree , obj_sample , obj_optCART , alpha0 )

% Maximum iteration for pruning
mI = nnodes( obj_tree ) + 1;

if nargin < 4
    alpha0 = Inf;
end

% Pre-set
alpha   = NaN( 1 , mI );
prTrees = repmat( tree , [ 1 , mI ] );
iOS     = NaN( 1 , mI );

% Set first values
prTree        = obj_tree;
alpha( 1 )    = 0;
prTrees( 1 )  = prTree;
p_l_ID        = findleaves( prTrees( 1 ) );
MSE_T         = obj_optCART.crit_ios( [ obj_tree.Node{ p_l_ID } ] , obj_sample , obj_optCART );
iOS( 1 )      = MSE_T;


for i = 1 : mI
    
    % Number of nodes found in the current tree
    nPT = nnodes( prTree );
    % If root node is found break the loop
    if nPT == 1
        break;
    end
    % Create candidate tree's -> find nodes to be chopped
    cand_leaves = findleaves( prTree );
    cand_IDs    = getparents( prTree , cand_leaves );
    cand_IDs    = unique( cand_IDs );
    nPT_j = numel( cand_IDs );
    % Create chopped trees for the most recent tree
    chTree    = repmat( tree , [ nPT_j , 1 ] );
    MSE_ch    = NaN( nPT_j , 1 );
    deltas    = MSE_ch;
    numLeaves = MSE_ch;
    for j = 1 : nPT_j
        % Chopping the recent tree
        chTree( j )    = chopTree( prTree, cand_IDs( j ) );
        % Get the in-sample criterion
        c_l_ID         = findleaves( chTree( j ) );
        numLeaves( j ) = numel( c_l_ID );
        MSE_ch( j )    = obj_optCART.crit_ios( [ chTree( j ).Node{ c_l_ID } ] , obj_sample , obj_optCART );
        % Get alpha candidates
        deltas( j ) = ( MSE_T - MSE_ch( j ) ) ./ ( numLeaves( j ) - nPT );
    end
    % Find the weakest node
    [ ~ , wD_ID ] = min( deltas );
    alpha_i = ( MSE_ch( j ) - iOS( 1 ) ) ./ numLeaves( j );
    if alpha_i > alpha0
        break;
    else
        % Set the new benchmark tree
        alpha( i + 1 ) = alpha_i;
        MSE_T   = MSE_ch( wD_ID );
        prTree  = chTree( wD_ID );
        prTrees( i + 1 ) = prTree;
        iOS( i + 1 )     = MSE_T;
    end
end

if i == mI
    error('Pruning did not converged!')
end

if nargin < 4
    alpha = alpha( 1 : i );
    prTrees = prTrees( 1 : i );
    iOS = iOS( 1 : i );
else
    alpha   = alpha0;
    prTrees = prTree;
    iOS     = MSE_T;
end

end