function [ y_hat , vars ] = infJack( predictions , inbag , used_trees )

num_tree = size( inbag , 2 );

if nargin < 3
    % Which trees to use
    used_trees = 1 : num_tree;
end

% Filter for the predictions
pred = predictions( : , used_trees );
% Check if bootstrapping is without replacement
no_replacement = all( max( inbag ) == 1 );

% Extract tree-wise predictions and variable counts from random forest

% number of tree
B = numel( used_trees );
% Number of observations
n = size( inbag , 1 );
% Average number of observations used
% s = sum( inbag( : ) ) / size( inbag , 2 );

% Average prediction for each observation
y_hat = mean( pred , 2 );
% Center the predictions
pred_centered = pred - y_hat;

% Create sparse matrix
N = sparse( inbag( : , used_trees ) );
% Average of usage of each observations
N_avg = mean( inbag , 2 );

if ( B ^ 2 > n * size( pred , 1 ) )
    C = N * pred_centered' - N_avg * sum( pred_centered , 2 )';
    raw_IJ = ( sum( C.^2 ) ./ B.^2 )';
else
    NTN = N' * N;
    NTNPT_T = pred_centered * NTN';
    T1 = sum( pred_centered .* NTNPT_T , 2 );
    
    RS = sum( pred_centered , 2 );
    NbarTN = N_avg' * N;
    T2 = RS .* ( NbarTN * pred_centered' )';
    
    T3 = sum( N_avg(:).^2 ).* RS.^2;
    raw_IJ = ( T1 - 2 .* T2 + T3 ) ./ B.^2;
end

% Apply Monte Carlo bias correction

N_var = mean( mean( N .^ 2 , 2 ) - mean( N , 2 ).^ 2 );
boot_var = sum( pred_centered .^ 2 , 2 ) ./ B;
bias_correction = n .* N_var .* boot_var ./ B;
vars = raw_IJ - bias_correction;

% Finite sample correction
if no_replacement
    variance_inflation = 1 ./ ( 1 - mean( inbag( : ) ) ) .^ 2;
    vars = variance_inflation .* vars;
end


end