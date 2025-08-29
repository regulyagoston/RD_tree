%% Creates the splitting values for Zi
%  Agoston Reguly
%


function [ bS_f , valid_split , spQ_i , split , nS ] = createSplits( Zi , ITTi , node , tr_tree , id_f , parentID , data , optTree )

% Sorted and unique values of the splitting variable
Z_u = unique( Zi ); 
% Number of potential split
nZ = numel( Z_u );
%% Get general splits created by the split functions:
% 1st type: 
%   If the potential splits are less or equal than the maximum number of split,
%   then take each value as a splitting value
if nZ <= optTree.numSplit
    split = Z_u;
    split( end ) = [];
else
%   Other specified splitting criterions which are using only the
%   feature and the number of splits
    split = optTree.split_fnc( Zi , optTree.numSplit );
    split( end ) = [];
end


%% Check for the validity of the split - enough observations in each split option
% Number of potential split
ns_s = numel( split );
% Logical for validity of the split
log_valid = false( ns_s , 1 );
% logical for isempty Wi
lWi = isempty( ITTi );
if ~lWi
    nITTi = ~ITTi;
end

estr_i = tr_tree.Node{ parentID }.est_struct;
log_bw = isfield( estr_i, 'h_j' );

for i = 1 : ns_s
    % Logicals for below and above the split value
    log_above = Zi > split( i );
    log_below = Zi <= split( i );
    if lWi
        % If no treatment
        nA = sum( log_above );
        nB = sum( log_below );
        % True if both below and above parts have enough observations
        log_valid( i ) = nA >= optTree.minObs && nB >= optTree.minObs;
    else
        
        if log_bw % check the effective number of observations with bandwidth
            X_i = data.X_tr( node.logID_tr );
            kernFnc = str2func( [ optTree.kernel , '2' ] );
            % Running variable and wieghts BELOW treshold and BELOW split point
            Xc_b = X_i( log_below & nITTi , : );
            Wc_b = feval( kernFnc , ( ( Xc_b - data.c ) ./ estr_i.h_j ) ) ./ estr_i.h_j;
            nC_b = sum( Wc_b ~= 0 );
            % Running variable and wieghts BELOW treshold and ABOVE split point
            Xc_a = X_i( log_above & nITTi , : );
            Wc_a = feval( kernFnc , ( ( Xc_a - data.c ) ./ estr_i.h_j ) ) ./ estr_i.h_j;
            nC_a = sum( Wc_a ~= 0 );
            % Running variable and wieghts ABOVE treshold and BELOW split point
            Xt_b = X_i( log_below & ITTi , : );
            Wt_b = feval( kernFnc , ( ( Xt_b - data.c ) ./ estr_i.h_j ) ) ./ estr_i.h_j;
            nT_b = sum( Wt_b ~= 0 );
            % Running variable and wieghts ABOVE treshold and ABOVE split point
            Xt_a = X_i( log_above & ITTi , : );
            Wt_a = feval( kernFnc , ( ( Xt_a - data.c ) ./ estr_i.h_j ) ) ./ estr_i.h_j;
            nT_a = sum( Wt_a ~= 0 );
        else
            % Number of control units below and above the split
            nC_b = sum( log_below & nITTi );
            nC_a = sum( log_above & nITTi );
            % Number of treatment units below and above the split
            nT_b = sum( log_below & ITTi );
            nT_a = sum( log_above & ITTi );
        end
        % True if both treatment and control has enough observations
        log_valid( i ) = nC_b >= optTree.minObs && nC_a >= optTree.minObs && ...
                         nT_b >= optTree.minObs && nT_a >= optTree.minObs;
    end
end

%% Create valid splits
split = split( log_valid );
nS = numel( split );


%% Compare candidates by each split using splitting criterion
if nS > 0 
    [ bS_f , spQ_i ] = get_best_split( node , tr_tree , split , nS , id_f , parentID , data , optTree );
    valid_split = true;
    split = bS_f( 1 ).splitVal( end );
else
    split                     = NaN;
    spQ_i                     = NaN;
    bS_f                      = node;
    bS_f.Failed               = 2;
    bS_f.feature_chck( id_f ) = false;
    valid_split               = false;
end

end


function [ bS_f , spQ_i ] = get_best_split( node , tr_tree , split , nS , id_f , parentID , data , optTree )

%% Candidate Nodes
c_nodes  = repmat( nodeProp , [ nS , 2 ] );
nSubSet  = [ node.subsetIdx , id_f ];
nLEq     = [ node.LEq , true ];
nGE      = [ node.LEq , false ];

%% Update all candidate
for i = 1 : nS
    % Create the nodes
    nSV  = [ node.splitVal , split( i ) ];
    c_nodes( i , 1 ) = nodeProp( nSubSet , nSV , nLEq , 0 , node.feature_chck );
    c_nodes( i , 2 ) = nodeProp( nSubSet , nSV , nGE  , 0 , node.feature_chck );
    % Update them -> using only the training sample
    c_nodes( i , 1 ) = updateNode( c_nodes( i , 1 ) , data , optTree , true );
    c_nodes( i , 2 ) = updateNode( c_nodes( i , 2 ) , data , optTree , true );
end

%% Find the smallest criterion value among the candidates
if nS > 1
    % Calculate the splitting criterion for each candidate
    spQ_i = NaN( nS , 1 );
    lfID = findleaves( tr_tree );
    cfID = lfID( lfID ~= parentID );
    for i = 1 : nS
        cand_tree = [ tr_tree.Node{ cfID } , c_nodes( i , : ) ];
        spQ_i( i ) = optTree.crit_split( cand_tree , data , optTree );
    end
    % Choose the smallest as best split
    [ ~ , minId ] = min( spQ_i );
    bS_f          = c_nodes( minId , : );
else
    % If only one candidate
    bS_f  = c_nodes;
    spQ_i = NaN;
end

end