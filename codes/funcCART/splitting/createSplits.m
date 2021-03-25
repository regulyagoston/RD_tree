%% Creates the splitting values for Zi
%  Agoston Reguly
%


function [ bS_f , valid_split , spQ_i , split , nS ] = createSplits( Zi , Wi , node , tr_tree , id_f , parentID , data , optTree )

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
    split = optTree.splitFnc( Zi , optTree.numSplit );
    split( end ) = [];
end


%% Check for the validity of the split - enough observations in each split option
% Number of potential split
ns_s = numel( split );
% Logical for validity of the split
log_valid = false( ns_s , 1 );
% logical for isempty Wi
lWi = isempty( Wi );
if ~lWi
    nWi = ~Wi;
end
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
        % Number of control units below and above the split
        nC_b = sum( log_below & nWi );
        nC_a = sum( log_above & nWi );
        % Number of treatment units below and above the split
        nT_b = sum( log_below & Wi );
        nT_a = sum( log_above & Wi );
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