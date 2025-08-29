%% Cross validate Random Forest parameters:
%   number of trees
%   number of variables to choose at splits
%   bandwidth parameter


function [ forest_fin, inbag_fin, ...
    num_tree, mtry, bw, oOB_crit, ...
    forest, inbag ] = cross_validate_forest( obj_sample, optTree, num_tree, mtry, bw )

if isempty( num_tree )
    num_tree = 500;
    nTree = 1;
else
    nTree = numel( num_tree );
end

if isempty( mtry )
    mtry = get_num_vars_forest( size(obj_sample.Z_tr,2), [] );
    nMtry = 1;
else
    nMtry = numel( mtry );
end

if isempty( bw ) && strcmp( optTree.model, 'local' )
    bw = bw_grid(obj_sample, obj_optCART);
    nBW = numel( bw );
elseif ~isempty( bw ) && strcmp( optTree.model, 'local' )
    nBW = numel( bw );
else 
    nBW = 1;
    bw = NaN;
end


oOB_crit = NaN( nTree, nMtry, nBW );
forest = repmat( tree , [ nTree, nMtry, nBW ] );
inbag = cell( [ nTree, nMtry, nBW ] );
optTree_ijk = optTree;
for i = 1:nTree
    for j = 1:nMtry
        optTree_ijk.mtry = mtry( j );
        for k = 1:nBW
            optTree_ijk.bw = bw( k );
            [ forest( i, j, k ) , oOB_crit( i, j, k ), inbag_aux ] = runForest( obj_sample , optTree_ijk , num_tree( i ) );
            inbag{ i, j, k } = inbag_aux;
        end
    end
end

% Find the minimum value and its linear index
[~, linearIdx] = min(oOB_crit(:));

% Convert linear index to subscript indices
[row, col, page] = ind2sub(size(oOB_crit), linearIdx);

forest_fin = forest( row, col, page );
inbag_fin = inbag{ row, col, page };





end