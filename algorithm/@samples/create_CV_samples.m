%% Create K-Fold subsamples which relates to samples



function [ train_sample , test_sample ] = create_CV_samples( obj_samples , obj_optCART , K )

%% Error check for the functions
if ~isa( obj_samples , 'samples' )
    error( 'samples:create_CV_samples:wrongInput:obj_samples',...
           'Invalid input for obj_samples, use samples object!' )
end
if ~isa( obj_optCART , 'optCART' )
    error( 'samples:create_CV_samples:wrongInput:obj_optCART',...
           'Invalid input for obj_optCART, use optCART object!' )
end
if K < 2 || mod( K , 1 ) ~= 0
    error( 'samples:create_CV_samples:wrongInput:K',...
           'Invalid number of cross-validation sample use a positive integer, larger than 1!' )
end
   
%% Create the partition indexes
% 1st get the partition indexes
kIdx_tr = createPartition( obj_samples.n_tr , K , obj_optCART.seed );
if obj_optCART.honest
    kIdx_est = createPartition( obj_samples.n_est , K , obj_optCART.seed + 1 );
end


%% Create the cross-validation samples
train_sample = [];
test_sample  = [];
for k = 1 : K
    % Set up the new inputs
    var_tr = {};
    var_te = {};
    
    % TRAINING: Choose the indexes for the K'th sample (all except k'th column)
    id_tr_tr_k = [ kIdx_tr( : , 1 : k - 1 ) , kIdx_tr( : , k + 1 : end ) ];
    id_tr_tr_k = unique( id_tr_tr_k( : ) , 'stable' );
    % TESTING: k'th column
    id_te_tr_k = kIdx_tr( : , k );
    % In case of honest method
    if obj_optCART.honest
        % TRAINING sample's estimation sample
        id_tr_est_k = [ kIdx_est( : , 1 : k - 1 ) , kIdx_est( : , k + 1 : end ) ];
        id_tr_est_k = unique( id_tr_est_k( : ) , 'stable' );
        % TESTING sample's estimation sample
        id_te_est_k = kIdx_est( : , k );
    end
    
    % Get variables Y and Z and training/estimation indexes
    if obj_optCART.honest
        %% Samples - tr samples
        Y_tr_tr_k  = obj_samples.Y_tr(  id_tr_tr_k  , : );
        Y_tr_est_k = obj_samples.Y_est( id_tr_est_k , : );
        n_tr_tr_k  = numel( Y_tr_tr_k );
        n_tr_est_k = numel( Y_tr_est_k );
        Z_tr_tr_k  = obj_samples.Z_tr(  id_tr_tr_k  , : );
        Z_tr_est_k = obj_samples.Z_est( id_tr_est_k , : );
        Y_tr = [ Y_tr_tr_k ; Y_tr_est_k ];
        Z_tr = [ Z_tr_tr_k ; Z_tr_est_k ];
        index4TrEst_tr = make_index4TrEst( n_tr_tr_k , n_tr_est_k );
        var_tr = [ var_tr , { 'index4TrEst' , index4TrEst_tr } ];
        
        % Samples - te samples
        Y_te_tr_k  = obj_samples.Y_tr(  id_te_tr_k  , : );
        Y_te_est_k = obj_samples.Y_est( id_te_est_k , : );
        n_te_tr_k  = numel( Y_te_tr_k );
        n_te_est_k = numel( Y_te_est_k );
        Z_te_tr_k  = obj_samples.Z_tr(  id_te_tr_k  , : );
        Z_te_est_k = obj_samples.Z_est( id_te_est_k , : );
        Y_te = [ Y_te_tr_k ; Y_te_est_k ];
        Z_te = [ Z_te_tr_k ; Z_te_est_k ];
        index4TrEst_te = make_index4TrEst( n_te_tr_k , n_te_est_k );
        var_te = [ var_te , { 'index4TrEst' , index4TrEst_te } ];
        
    else
        Y_tr = obj_samples.Y_tr(  id_tr_tr_k , : );
        Z_tr = obj_samples.Z_tr(  id_tr_tr_k , : );
        Y_te = obj_samples.Y_tr(  id_te_tr_k , : );
        Z_te = obj_samples.Z_tr(  id_te_tr_k , : );
    end
    
    % Add treatment dummy if needed
    if ( strcmp( obj_optCART.type , 'CATE' ) && ~strcmp( obj_optCART.criteria , 'RDD' ) ) ...
            || strcmp( obj_optCART.type , 'CLATE' )
        if obj_optCART.honest
            W_tr = [ obj_samples.W_tr(  id_tr_tr_k  , : ) ;...
                     obj_samples.W_est( id_tr_est_k , : ) ];
            W_te = [ obj_samples.W_tr(  id_te_tr_k  , : ) ;...
                     obj_samples.W_est( id_te_est_k , : ) ];
            W_tr = [ obj_samples.W_tr(  id_tr_tr_k  , : ) ;...
                     obj_samples.W_est( id_tr_est_k , : ) ];
            W_te = [ obj_samples.W_tr(  id_te_tr_k  , : ) ;...
                     obj_samples.W_est( id_te_est_k , : ) ];
        else
            W_tr = obj_samples.W_tr( id_tr_tr_k , : );
            W_te = obj_samples.W_tr( id_te_tr_k , : );
        end
        % Convert to logicals
        W_tr = W_tr == 1;
        var_tr = [ var_tr , { 'W' , W_tr } ];
        W_te = W_te == 1;
        var_te = [ var_te , { 'W' , W_te } ];
    end
    
    % Extract the running variable/explanatory variables
    if strcmp( obj_optCART.model , 'linear' ) || strcmp( obj_optCART.criteria , 'RDD' )
        % If linear or RDD is used add X and clusters
        sX = obj_samples.nCov - obj_optCART.orderPolinomial;
        if sX ~= 0
            % Select the given covariates
            idX = [ 1 , obj_optCART.orderPolinomial + 1 : obj_samples.nCov ];
        else
            idX = 1;
        end
        if obj_optCART.honest
            X_tr = [ obj_samples.X_tr(  id_tr_tr_k , idX ) ; ...
                     obj_samples.X_est( id_tr_est_k , idX ) ];
            X_te = [ obj_samples.X_tr(  id_te_tr_k , idX ) ; ...
                     obj_samples.X_est( id_te_est_k , idX ) ];
            if isprop( obj_samples , 'cl_tr' )
                cl_tr = [ obj_samples.cl_tr(  id_tr_tr_k , : ) ; ...
                          obj_samples.cl_est( id_tr_est_k , : ) ];
                cl_te = [ obj_samples.cl_tr(  id_te_tr_k , : ) ; ...
                          obj_samples.cl_est( id_te_est_k , : ) ];
            end
        else
            X_tr = obj_samples.X_tr(  id_tr_tr_k , idX );
            X_te = obj_samples.X_tr(  id_te_tr_k , idX );
            if isprop( obj_samples , 'cl_tr' )
                cl_tr = obj_samples.cl_tr(  id_tr_tr_k , : );
                cl_te = obj_samples.cl_tr(  id_te_tr_k , : );
            end
        end
        var_tr = [ var_tr , { 'X' , X_tr } ];
        var_te = [ var_te , { 'X' , X_te } ];
        if strcmp( obj_optCART.model , 'linear' ) && strcmp( obj_optCART.varEstimator , 'clustered' )
            var_tr = [ var_tr , { 'cluster' , cl_tr } ];
            var_te = [ var_te , { 'cluster' , cl_te } ];
        end
    end
    if strcmp( obj_optCART.criteria , 'RDD' )
        var_tr = [ var_tr , { 'cutoff' , obj_samples.c } ];
        var_te = [ var_te , { 'cutoff' , obj_samples.c } ];
    end
    
    % Create the cross-validation sample
    aux_sample_tr = samples( obj_optCART , Y_tr , Z_tr , var_tr{ : } );
    train_sample  = [ train_sample ; aux_sample_tr ];
    aux_sample_te = samples( obj_optCART , Y_te , Z_te , var_te{ : } );
    test_sample  = [ test_sample ; aux_sample_te ];
end
            
end


%% Create the indexes for the training and estimation sample
function index4TrEst = make_index4TrEst( n_tr , n_est )

nM = max( n_tr , n_est );
if nM == n_tr && nM == n_est
    index4TrEst = [ 1 : n_tr ; n_tr + 1 : n_tr + n_est ]';
elseif nM > n_est
    index4TrEst = [ 1 : n_tr ; [ n_tr + 1 : n_tr + n_est , ...
        repmat( n_tr + n_est , [ 1 ,  2 * nM - ( n_est + n_tr ) ] ) ] ]';
else
    index4TrEst = [ [ 1 : n_tr , ...
        repmat( n_tr , [ 1 , nM - n_tr ] ) ] ; ...
        n_tr + 1 : n_tr + n_est ]';
end

end