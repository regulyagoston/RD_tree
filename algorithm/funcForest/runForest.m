function [ forest , oOS , inbag ] = runForest( sample , optTree , num_tree )
%RUNFOREST Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    num_tree = 500;
end

if isprop( sample , 'X_tr' )
    logX = true;
else
    logX = false;
end
if isprop( sample , 'W_tr' )
    logW = true;
else
    logW = false;
end
    
%% Create bootstrap samples
% Get the observations to use and save to inbag the frequencies
if optTree.bootstrap_w_replacement
    bs_Idx = createBootstrapSampleIdx( sample.Y_tr , sample.n_tr , num_tree , optTree.seed );
else
    no_bootstrap_sample = floor( optTree.sample_fraction .* sample.n_tr );
    bs_Idx = NaN( no_bootstrap_sample , num_tree );
    for i = 1 : num_tree
        rng( optTree.seed + i )
        bs_Idx( : , i ) = randperm( sample.n_tr , no_bootstrap_sample );
    end
end
inbag  = zeros( sample.n_tr , num_tree );

% for n = 1 : sample.n_tr
%     inbag( n , : ) = sum( bs_Idx == n );
% end
for c = 1:num_tree
    inbag(:, c) = accumarray(bs_Idx(:, c), 1, [sample.n_tr, 1]);
end
%% Run for each tree
gr_tree  = repmat( tree , [ num_tree , 1 ] );
forest   = repmat( tree , [ num_tree , 1 ] );
%poi_ij   = NaN( sample.n_tr , num_tree );
oOS      = NaN( num_tree, 1 );
 
optTree_i = copyWithChange(optTree,'forest', false);

all_idx = 1:sample.n_tr;
for i = 1 : num_tree
    % select features and observations
    Yi_tr = sample.Y_tr( bs_Idx( : , i ) , : );
    Zi_tr = sample.Z_tr( bs_Idx( : , i ) , : );
    % Select out-of-bag observations
    OoB_log = setdiff( all_idx, bs_Idx(:,i) );
    Yi_est = sample.Y_tr( OoB_log , : );
    Zi_est = sample.Z_tr( OoB_log , : );
    N_tr = numel( Yi_tr );
    N_est = numel( Yi_est );
    % Join for sample object
    Yi = [Yi_tr; Yi_est ];
    Zi = [Zi_tr; Zi_est ];   
    %Ni = numel( Yi );
    if logX
        Xi_tr = sample.X_tr( bs_Idx( : , i ) , : );
        Xi_est = sample.X_tr( OoB_log , : );
        Xi = [Xi_tr; Xi_est ];
    end
    if logW
        Wi_tr = sample.W_tr( bs_Idx( : , i ) , : );
        Wi_est = sample.W_tr( OoB_log , : );
        Wi = [Wi_tr; Wi_est ];
    end    
    
    idxTrEst = [ 1 : N_tr; ...
        [ ( N_tr + 1 ) : ( N_tr + N_est ), repmat( ( N_tr + 1 ), [1,N_tr - N_est]) ] ]';

    % Set the sample which is dependent on the method
    if logX & ~logW
        obj_sample_i = samples( optTree_i , Yi  , Zi , 'index4TrEst' , idxTrEst ,...
                                      'X' , Xi , 'cutoff' , sample.c );
    elseif logW & ~logX
        obj_sample_i = samples( optTree_i , Yi  , Zi , 'index4TrEst' , idxTrEst, ...
                                       'W' , Wi );
    elseif logW & logX
        obj_sample_i = samples( optTree_i , Yi  , Zi , 'index4TrEst' , idxTrEst, ...
                                       'W' , Wi, 'X', Xi, 'cutoff' , sample.c );
    else
        obj_sample_i = samples( optTree_i , Yi  , Zi , 'index4TrEst' , idxTrEst );
    end

    % Grow i'th tree
    gr_tree( i ) = growRDDtree( obj_sample_i , optTree_i );
    
    % Estimate i'th tree on the estimation sample
    forest( i ) = estimate_tree_est( gr_tree( i ) , obj_sample_i , optTree_i );
    
    % Predict point of interest for each observations
    % poi_ij( : , i ) = predict( forest( i ) , sample.Z_tr ); 

    % Calculate out of bag criterion
   p_l_ID = findleaves( forest( i ) );
   oOS( i )  = optTree_i.crit_oos( [ forest( i ).Node{ p_l_ID } ] , obj_sample_i , optTree_i );
                     
end

% Make prediction from poi_ij and get the variances
% [ predicted , poi_var ] = forestInfJack( poi_ij , inbag , true , 1 : num_tree );


end

