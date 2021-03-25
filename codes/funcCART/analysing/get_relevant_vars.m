%% Get relevant features from a tree



function [ var_id , rel_var_c , split_val ] = get_relevant_vars( obj_tree )

lID = findleaves( obj_tree );
nL = numel( lID );
rel_vars = [];
split_val = cell( nL );
rel_var_c = cell( nL );

for i = 1 : nL
    rel_var_c{ i } = obj_tree.Node{ lID( i ) }.subsetIdx;
    split_val{ i } = obj_tree.Node{ lID( i ) }.splitVal;   
    rel_vars  = [ rel_vars , rel_var_c{ i } ]; 
end

% Get unique variable indexes
var_id = unique( rel_vars );


end