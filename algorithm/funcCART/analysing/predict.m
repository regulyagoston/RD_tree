%% Predict the outcome along the feature values
% Input:
%   tree - estimated tree structure
%   Z    - corresponding features Z to be evaluated
%  Optional: 'struct_val':
%       char -> value to predict:
%
% Output:
%   pred - predicted outcome for the values of Z
%
% Last updated: 29/04/2021


function pred = predict( tree , Z , struct_val )

if nargin < 3
    % Make an inteligent guess for the predicted outcome based on existing
    % fields
    if isfield( tree.Node{ 1 }.est_struct , 'tau_i' )
        struct_val = 'tau_i';
    elseif isfield( tree.Node{ 1 }.est_struct , 'Y_hat_i' )
        struct_val = 'Y_hat_i';
    else
       error('Not implemented the automatic prediction for such tree!') 
    end
end

if ~isfield( tree.Node{ 1 }.est_struct , struct_val )
    error(horzcat( struct_val , ' is not an implemented prediction value for the estimated tree!' ) )
end

% Set values
N = size( Z , 1 );
pred = NaN( N , 1 );

% Find leaves
lID = findleaves( tree );
nL = numel( lID );

if any( strcmpi( struct_val( end-1 : end ) , '_i' ) )
    log_i = true;
else
    log_i = false;
end

% Get the treatment values
for i = 1 : nL
    % Select the leaf
    node_j = tree.Node{ lID( i ) };
    % Get logicals for the covariates
    leaf_id_j = get_logIdx( node_j , Z );
    % Save the predicted outcome
    aux = node_j.est_struct.( struct_val );
    if log_i
        pred( leaf_id_j )    = aux( 1 );
    else
        pred( leaf_id_j )    = aux;
    end
end

end