%% Cross validate to get optimal pruning/complexity parameter


function [ opt_beta , s_oOS , beta , oOS ] = cross_validate( obj_tree , obj_sample , obj_optCART )


% 1st: Get the raw complexity parameters
[ prTrees , iOS , alpha ] = pruning( obj_tree , obj_sample , obj_optCART );

% 2nd: calculate 'typical values' instead of raw values
nM = numel( alpha );
beta = zeros( 1 , nM );
if any( alpha < 0 )
    %warning('There are negative candidate complexity parameters!')
    signA = sign( alpha );
else
    signA = ones( 1 , nM );
end
for i = 2 : nM - 1
    beta( i ) =  sqrt( abs( alpha( i ) .* alpha( i + 1 ) ) );
end
beta( end ) = Inf;
beta = beta .* signA;
if numel( unique( beta ) ) ~= numel( beta )
    error('cross_validate','Some problems with the complexity parameters! There are repetitions of the same value!')
end

%3rd: create cross-validation samples
[ train_sample , test_sample ] = create_CV_samples( obj_sample , obj_optCART , obj_optCART.numKfold );
oOS = NaN( obj_optCART.numKfold , nM );
% Prune with betas and predict on test sample
prTrees_k   = repmat( tree , [ obj_optCART.numKfold , nM ] );
predTrees_k = repmat( tree , [ obj_optCART.numKfold , nM ] );
try
    if ~obj_optCART.paralell
       error('This is only a shortcut if paralell processing is not allowed!') 
    end
    parfor k = 1 : obj_optCART.numKfold
        % Fit a full model on the cv data set
        obj_tree_k = growRDDtree( train_sample( k ) , obj_optCART );
    %     % Old codes:
    %     %obj_tree_k = clearTree( obj_tree );
    %     % Update with the k'th training sample
    %     %obj_tree_k = updateTree( obj_tree_k , [] , train_sample( k ) , obj_optCART , true );
        % For each candidate complexity parameter
        for i = 1 : nM
            % Prune trees based on training sample
            prTrees_k( k , i )   = pruning( obj_tree_k , train_sample( k ) , obj_optCART , beta( i ) );
            % Use the pruned trees to predict or estimate on test sample
            predTrees_k( k , i ) = clearTree( prTrees_k( k , i ) );
            if any( strcmp( obj_optCART.type , 'prediction' ) )
                predTrees_k( k , i ) = predictTree( predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k ) , obj_optCART );
            elseif any( strcmp( obj_optCART.type , {'CATE','CLATE'} ) )
                predTrees_k( k , i ) = updateTree(  predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k )   , obj_optCART , true );
            else
                error('cross_validate:wrongInputType',['There are no such type implemeted as: ' obj_optCART.type ] )
            end
            % Calcilate out of sample criterion on test sample
            p_l_ID        = findleaves( predTrees_k( k , i ) );
            oOS( k , i )  = obj_optCART.crit_oos( [ predTrees_k( k , i ).Node{ p_l_ID } ] , test_sample( k ) , obj_optCART );
        end
    end
catch
    for k = 1 : obj_optCART.numKfold
        % Fit a full model on the cv data set
        obj_tree_k = growRDDtree( train_sample( k ) , obj_optCART );
    %     % Old codes:
    %     %obj_tree_k = clearTree( obj_tree );
    %     % Update with the k'th training sample
    %     %obj_tree_k = updateTree( obj_tree_k , [] , train_sample( k ) , obj_optCART , true );
        % For each candidate complexity parameter
        for i = 1 : nM
            % Prune trees based on training sample
            prTrees_k( k , i )   = pruning( obj_tree_k , train_sample( k ) , obj_optCART , beta( i ) );
            % Use the pruned trees to predict or estimate on test sample
            predTrees_k( k , i ) = clearTree( prTrees_k( k , i ) );
            if any( strcmp( obj_optCART.type , 'prediction' ) )
                predTrees_k( k , i ) = predictTree( predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k ) , obj_optCART );
            elseif any( strcmp( obj_optCART.type , {'CATE','CLATE'} ) )
                predTrees_k( k , i ) = updateTree(  predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k )   , obj_optCART , true );
            else
                error('cross_validate:wrongInputType',['There are no such type implemeted as: ' obj_optCART.type ] )
            end
            % Calcilate out of sample criterion on test sample
            p_l_ID        = findleaves( predTrees_k( k , i ) );
            oOS( k , i )  = obj_optCART.crit_oos( [ predTrees_k( k , i ).Node{ p_l_ID } ] , test_sample( k ) , obj_optCART );
        end
    end
end

% Get the sum of out-of-sample values
s_oOS = sum( oOS );
% Get the minimum value
[ mVal , mID ] = min( s_oOS );
% If 1SE rule
if obj_optCART.CV1SE
    u_range  = mVal + ( std( oOS( : , mID ) ) ./ obj_optCART.numKfold );
    val_id   = s_oOS <= u_range;
    aux_beta = beta( val_id );
    opt_beta = aux_beta( end );
else
    opt_beta = beta( mID );
end



end