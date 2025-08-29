%% Find the best split for a given feature
% It works with the training sample only


function [ bS_f , valid_split , cand_IOS ] = best_split_by_feature( tr_tree , parentID , id_f , data , optTree )

% Initialization
valid_split = true;
node        = tr_tree.Node{ parentID };

%% Get the selected feature and the possible treatment variable
Zi = data.Z_tr( node.logID_tr , id_f );
Zu = unique( Zi );
nZ = numel( Zu );

Yi = data.Y_tr( node.logID_tr , : );
if data.log_W_tr
    % Get units who are treated within that leaf
    Wi = data.W_tr( node.logID_tr , : );
else
    Wi = [];
end
if data.log_ITT_tr
    % Get units who are treated within that leaf
    ITTi = data.ITT_tr( node.logID_tr , : );
else
    ITTi = [];
end
if data.log_X_tr
    % Get controls within that leaf
    Xi = data.X_tr( node.logID_tr , : );
else
    Xi = [];
end
if data.log_cl_tr
    cl_i = data.cl_tr( node.logID_tr , : ); 
else
    cl_i = NaN;
end
if optTree.honest
    % Get outcome variables from estimation sample to control for enough
    % observations
    W_est = data.W_est( node.logID_est , : );
    Z_est = data.Z_est( node.logID_est , id_f );
    ITT_est = data.ITT_est( node.logID_est , : );
end

%% Find best split for the given feature
%   - Heart of the algorithm
%   - Different criterion uses different splitting rules

%%% Different implemented (fast) rutines:
bf = false;
warnState = warning('off', 'all');
% 1st -> rpart prediction, using means and MSE leaf-by-leaf
if strcmp( optTree.splitType , 'rpart' ) && strcmp( optTree.type , 'prediction' ) ...
   && strcmp( optTree.model , 'mean' )   && strcmp( optTree.criteria , 'MSE' )

   split = split_rpart_classical( Yi , Zi , optTree );
   
% 2nd -> Athey-Imbens honest criterion with unconfoundedness
elseif nZ > optTree.numSplit && ... 
       any( strcmp( optTree.splitType , {'causal-tree-AI','causal-tree-RDD','causal-tree-npRDD-uni','causal-tree-npRDD-lbl'} ) )
    % Sort the variables according to the splitting value
    [ Zs , sort_ID ]       = sort( Zi );
    [ Z_e  , sort_ID_est ] = sort( Z_est );
    
    
    if strcmp( optTree.splitType , 'causal-tree-AI' )
        split = split_CT_AI( Yi( sort_ID , : ) , Zs , Wi( sort_ID , : ) , Z_e , W_est( sort_ID_est , : ) , data.train_to_est_ratio , optTree );        
    elseif strcmp( optTree.splitType , 'causal-tree-RDD' )
        % Adjust the clusters as well if given
        if ~isnan( cl_i )
            cl_is = cl_i( sort_ID );
        else
            cl_is = NaN( numel( sort_ID ) , 1 );
        end
        split = split_CT_RDD( Yi( sort_ID , : ) , Wi( sort_ID , : ) , Zs , Xi( sort_ID , : ) , ITTi( sort_ID , : ) , cl_is , Z_e , ITT_est( sort_ID_est , : ) , optTree );  
    elseif strcmp( optTree.splitType , 'causal-tree-npRDD-uni' )
        % Adjust the clusters as well if given
        if ~isnan( cl_i )
            cl_is = cl_i( sort_ID );
        else
            cl_is = NaN( numel( sort_ID ) , 1 );
        end
        Xest = data.X_est( node.logID_est , : );
        split = split_CT_npRDD_uni( Yi( sort_ID , : ) , Wi( sort_ID , : ) , Zs , Xi( sort_ID , : ) , ITTi( sort_ID , : ) , cl_is , Z_e , ITT_est( sort_ID_est , : ) , Xest( sort_ID_est, : ) , data.c, data.n_tr, data.n_est, optTree );  
    elseif strcmp( optTree.splitType , 'causal-tree-npRDD-lbl' )
        % Adjust the clusters as well if given
        if ~isnan( cl_i )
            cl_is = cl_i( sort_ID );
        else
            cl_is = NaN( numel( sort_ID ) , 1 );
        end
        Xest = data.X_est( node.logID_est , : );
        split = split_CT_npRDD_lbl( Yi( sort_ID , : ) , Wi( sort_ID , : ) , Zs , Xi( sort_ID , : ) , ITTi( sort_ID , : ) , cl_is , Z_e , ITT_est( sort_ID_est , : ) , Xest( sort_ID_est, : ) , data.c, data.n_tr, data.n_est, optTree );  
    else
        error('No such splitting criterion implemented!')
    end
    
else
    % Brute force splitting criterion
    bf = true;
    [ bS_f , valid_split , ~ , split ] = createSplits( Zi , ITTi, node , tr_tree , id_f , parentID , data , optTree );
    
end
warning(warnState);                 % Restore original warning state
% Create the best split
if ~isnan( split )% && ~bf
    bS_f  = repmat( nodeProp , [ 1 , 2 ] );
    bS_f( 1 , 1 ) = nodeProp( [ node.subsetIdx , id_f ] , ...
                              [ node.splitVal , split ] , ...
                              [ node.LEq , true ] , 0 , node.feature_chck );
    bS_f( 1 , 2 ) = nodeProp( [ node.subsetIdx , id_f ] , ...
                              [ node.splitVal , split ] , ...
                              [ node.LEq , false ] , 0 , node.feature_chck );
    % Update
    bS_f( 1 , 1 ) = updateNode( bS_f( 1 , 1 ) , data , optTree , true );
    bS_f( 1 , 2 ) = updateNode( bS_f( 1 , 2 ) , data , optTree , true );
elseif isnan( split ) || ~bf
    % If there is no valid splitting value remained it is not a valid feature!
    bS_f                      = node;
    bS_f.Failed               = 2;
    bS_f.feature_chck( id_f ) = false;
    valid_split               = false;
    cand_IOS                  = Inf;
    return;
end

%% Calculate the IOS criterion
lfID      = findleaves( tr_tree );
cfID      = lfID( lfID ~= parentID );
cand_tree = [ tr_tree.Node{ cfID } , bS_f ];
cand_IOS  = optTree.crit_ios( cand_tree , data , optTree );

%% Check whether IOS is better than root tree
orig_tree = [ tr_tree.Node{ lfID } ];
orig_IOS = optTree.crit_ios( orig_tree , data , optTree );
% If the original has bettes IOS value than invalid node
%penalty_ios = ( numel( cfID ) + 2 ) .* optTree.cp;
rel_ios = ( orig_IOS - cand_IOS ) ./ abs( orig_IOS ); %cand_IOS ./ ios_0;
if rel_ios <= optTree.cp || isnan( cand_IOS )%( 1 - rel_ios ) <= penalty_ios
    bS_f                      = node;
    bS_f.Failed               = 1;
    bS_f.feature_chck( id_f ) = false;
    % Decide - it might be useful in later split!
    valid_split               = false;
end

end

