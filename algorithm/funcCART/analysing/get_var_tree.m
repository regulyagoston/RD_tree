%% Function to estimate CATE variance
% Collects the required values from the tree leaves
%


function var_cate = get_var_tree( obj_tree , obj_sample , tr )


leaf_id = findleaves( obj_tree );
leaf_nodes = [ obj_tree.Node{ leaf_id } ];
nL = numel( leaf_id );

% Output matrix
tau   = NaN( nL , 1 );
n_ell = NaN( nL , 1 );


% Collect inputs
for i = 1 : nL
    if tr
        n_ell( i ) = leaf_nodes( i ).n_j_tr;
    else
        n_ell( i ) = leaf_nodes( i ).n_j_est;
    end
    tau( i ) = leaf_nodes( i ).est_struct.tau_j;
end

% leaf probabilities
if tr
    N = obj_sample.n_tr;
else
    N = obj_sample.n_est;
end

leaf_prob = n_ell ./ N;

% Get the variance by Var decoposition
var_cate = sum( leaf_prob .* tau .^ 2 ) - ( sum( leaf_prob .* tau ) .^ 2 );



end
    
