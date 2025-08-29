
function [ y_hat , var_hat ] = forestInfJack( predictions , inbag , calibrate , used_trees )

if nargin < 3
    % Which trees to use
    used_trees = 1 : size( inbag , 2 );
    % Calibrate or not
    calibrate = true;
elseif nargin < 4
    % Which trees to use
    used_trees = 1 : size( inbag , 2 );
end
% number of tree
B = numel( used_trees );

[ y_hat , vars ] = infJack( predictions , inbag , used_trees );

if ( numel ( vars ) <= 20 && calibrate )
    calibrate = false;
    warning('No calibration with n <= 20');
end

% If appropriate, calibrate variance estimates; this step in particular
%   ensures that all variance estimates wil be positive.

if calibrate
    % Compute variance estimates using half the trees
    calibration_ratio = 2;
    n_sample = ceil( B ./ calibration_ratio );
    rng( 777 );
    use_tree_Idx = randperm( numel( used_trees ) , n_sample );
    [ ~ , var_ss ] = infJack( predictions( : , used_trees ) , inbag , use_tree_Idx );
    
    % Use this second set of variance estimates to estimate scale of Monte Carlo noise
    sigma2_ss = mean( ( var_ss -  vars ) .^2 );
    delta = n_sample / B;
    sigma2 = ( delta .^ 2 + ( 1 - delta ) .^ 2 ) ./ ( 2 .* ( 1 - delta ) .^2 ) .* sigma2_ss;
    
    % Use Monte Carlo noise scale estimate for empirical Bayes calibration
    vars_calibrated = calibrateEB( vars , sigma2 );
    var_hat = vars_calibrated;
else
    var_hat = vars;
end
    
end






