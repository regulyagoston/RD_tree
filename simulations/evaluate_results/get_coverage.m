%% Get coverage ratio for the simulations



function [ coverage, bias ] = get_coverage( tree , Z , trueCATE , ci_val)

% Set values
N      = size( Z , 1 );
tau    = NaN( N , 1 );
tau_se = NaN( N , 1 );
ate    = NaN( N , 1 );

% Find leaves
lID = findleaves( tree );
nL = numel( lID );

% Get the treatment values
for i = 1 : nL
    % Select the leaf
    node_j = tree.Node{ lID( i ) };
    % Get logicals for the covariates
    leaf_id_j = get_logIdx( node_j , Z );
    % Set the treatment and se values
    tau( leaf_id_j )    = node_j.est_struct.tau_j;
    tau_se( leaf_id_j ) = node_j.est_struct.tau_se_j;
    % Calculate the average treatment effect based on trueCATE
    ate( leaf_id_j ) = mean( trueCATE( leaf_id_j ) );
end


%% Calculate the hit ratio:
% Hit ratio for confidence intervals
lratio = ( 1 - ci_val ) / 2;
CI_xs = norminv( [ lratio , ci_val + lratio ] , 0 , 1 );

hitCI = tau >= ( ate + CI_xs( 1 ) .* tau_se ) & ...
        tau <= ( ate + CI_xs( 2 ) .* tau_se );

% Get the coverage ratio
coverage = mean( hitCI );
bias = mean( ate - tau );

end

