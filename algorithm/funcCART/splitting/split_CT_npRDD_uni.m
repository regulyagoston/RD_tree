%% CT- nonparametric splitting criterion with unified bw
%
%

function split  = split_CT_npRDD_uni( Y , W , Z , X , ITT , cl , Z_est , ITT_est , X_est, c, N_tr, N_est, optTree )


%% Initialization
minObs  = optTree.minObs;
% Whether it is fuzzy design or not
fuzzy = strcmp( optTree.type , 'CLATE' );
uITT = ~ITT;
uITT_est    = ~ITT_est;
n_j = numel( Y );
n_t = sum( ITT );
n_c = sum( uITT );
Y_t  = Y(  ITT );
Y_c  = Y( uITT );
X_t = X(  ITT , : );
X_c = X( uITT , : );
cl_t = cl(  ITT );
cl_c = cl( uITT );

if fuzzy
    W_t = W(  ITT );
    W_c = W( uITT );
else
    W = NaN( n_j, 1);
end


%%%%%%%%%%
%% Check for each splitting candidate
%%%%%%%%%%
% To be efficient and fast only check those which satisfies:
%   there are enough observation in each split (lower and upper parts)
%   there are enough observation in each bucket (number of treated and control)
%
% I use Sherman-Morrison formula with recursive updating:
%   faster than lscov even if it needs to compute for each observation

%% First check for enough observations in the 'tails'
kernFnc = str2func( [ optTree.kernel , '2' ] );
w_tr = feval( kernFnc , ( ( X_t - c ) ./ optTree.bw ) ) ./ optTree.bw;
w_co = feval( kernFnc , ( ( X_c - c ) ./ optTree.bw ) ) ./ optTree.bw;
w_est_tr = feval( kernFnc , ( ( X_est(ITT_est) - c ) ./ optTree.bw ) ) ./ optTree.bw;
w_est_co = feval( kernFnc , ( ( X_est(uITT_est) - c ) ./ optTree.bw ) ) ./ optTree.bw;

w_tr_log = w_tr > 0;
n_tr_eff = sum( w_tr_log );
w_co_log = w_co > 0;
n_co_eff = sum( w_co_log );
w_est_tr_log = w_est_tr > 0;
n_est_tr_eff = sum( w_est_tr_log );
w_est_co_log = w_est_co > 0;
n_est_co_eff = sum( w_est_co_log );


% Number of effective units
left_tr     = cumsum(  w_tr_log );
left_co    = cumsum( w_co_log );
left_tr_e   = cumsum(  w_est_tr_log );
left_co_e  = cumsum( w_est_co_log );
right_tr    = n_tr_eff - left_tr;
right_co   = n_co_eff - left_co;
right_tr_e  = n_est_tr_eff - left_tr_e;
right_co_e = n_est_co_eff - left_co_e;

% Valid number of effective observations
id_valid_left_tr = find( left_tr == minObs , 1 , 'first' );
id_valid_right_tr = find( right_tr == minObs , 1 , 'last' );
id_valid_left_co = find( left_co == minObs , 1 , 'first' );
id_valid_right_co = find( right_co == minObs , 1 , 'last' );
if isempty( id_valid_left_tr ) || isempty( id_valid_right_tr ) || ...
    isempty( id_valid_left_co ) || isempty( id_valid_right_co )
    split = NaN;
    return;
end
cs_ITT = cumsum(ITT);
cs_uITT = cumsum(uITT);
logNumObs = cs_ITT >= id_valid_left_tr & ...
            cs_ITT <= id_valid_right_tr & ...
            cs_uITT >= id_valid_left_co & ...
            cs_uITT <= id_valid_right_co;

id_valid_left_tr_est = find( left_tr_e == minObs , 1 , 'first' );
id_valid_right_tr_est = find( right_tr_e == minObs , 1 , 'last' );
id_valid_left_co_est = find( left_co_e == minObs , 1 , 'first' );
id_valid_right_co_est = find( right_co_e == minObs , 1 , 'last' );

if isempty( id_valid_left_tr_est ) || isempty( id_valid_right_tr_est ) || ...
    isempty( id_valid_left_co_est ) || isempty( id_valid_right_co_est )
    split = NaN;
    return;
end

cs_ITT_est = cumsum(ITT_est);
cs_uITT_est = cumsum(uITT_est);
logNumObs_e = cs_ITT_est >= id_valid_left_tr_est & ...
            cs_ITT_est <= id_valid_right_tr_est & ...
            cs_uITT_est >= id_valid_left_co_est & ...
            cs_uITT_est <= id_valid_right_co_est;


Z_e_valid = Z_est( logNumObs_e );

if ~isempty( Z_e_valid )
    logNumObs = logNumObs & ...
            Z >= Z_e_valid( 1 ) & ...
            Z <= Z_e_valid( end );
else
    split = NaN;
    return;
end
               
if ~any( logNumObs & ITT  ) || ~any( logNumObs & uITT  )
    split = NaN;
    return;
end

% Get the first and last valid observations to iterate through
aux       = Z( logNumObs );
lNumObsZ  = Z >= aux( 1 ) & Z <= aux( end );
v_obs_tr  = Z( lNumObsZ &  ITT );
v_obs_co  = Z( lNumObsZ & uITT );
% Intent-to-treat values and controls
Ztr = Z( ITT , 1 );
Zco = Z( uITT , 1 );
% where does valid treated and control observations would start
f_obs_tr_aux  = max( find( v_obs_tr( 1 )   == Ztr , 1 , 'last' ) , minObs );
l_obs_tr_aux  = min( find( v_obs_tr( end ) == Ztr , 1 , 'first' ) , n_t - minObs );
f_obs_co_aux  = max( find( v_obs_co( 1 )   == Zco , 1 , 'last' ) , minObs );
l_obs_co_aux  = min( find( v_obs_co( end ) == Zco , 1 , 'first' ) , n_c - minObs );
% Get the indexes for overall observations
f_obs_all = find( Z == max( Ztr( f_obs_tr_aux ) , Zco( f_obs_co_aux ) ) , 1 , 'last' );
l_obs_all = find( Z == min( Ztr( l_obs_tr_aux ) , Zco( l_obs_co_aux ) ) , 1 , 'first' );
% Adjust the first and last observations for both tr and co
f_obs_tr = find( Ztr <= Z( f_obs_all ) , 1 , 'last' );
f_obs_co = find( Zco <= Z( f_obs_all ) , 1 , 'last' );
l_obs_tr = find( Ztr >= Z( l_obs_all ) , 1 , 'first' );
l_obs_co = find( Zco >= Z( l_obs_all ) , 1 , 'first' );

if any( f_obs_all >= l_obs_all || f_obs_tr >= l_obs_tr || f_obs_co >= l_obs_co )
    split = NaN;
    return;
end

% Let compute the parameter values for all the possible split values
n_v = l_obs_all - f_obs_all + 1;


% But we have to deal with the same values in Z:
% Number of iterations: unique values in valid Z
% All observations
Zuq = unique( Z( f_obs_all : l_obs_all ) );
Zuq = Zuq(2:end);

% If there are no unique variables, return
if isempty( Zuq ) 
    split = NaN;
    return;
end

% Check for number of splits
if ( numel( Zuq ) > optTree.numSplit )
    % Previously using percentiles -- not the best as different pattern may
    % be independent from occurance (or even more extreme values e.g.)
    % % Get frequencies
    % freq_Z = histcounts( Z( f_obs_all : l_obs_all ), Zuq );
    % cum_prob = cumsum( freq_Z ) ./ numel( Z( f_obs_all : l_obs_all ) );
    % % Requested percentages
    % th_pctg = cumsum( ones( optTree.numSplit, 1 ) ) ./ optTree.numSplit;
    % Zuq_iter = NaN( optTree.numSplit, 1 );
    % for i = 1:optTree.numSplit
    %     [~,idx] = min( abs( cum_prob - th_pctg( i ) ) );
    %     Zuq_iter( i ) = Zuq( idx );
    % end
    % In case of closest pctgs remains same

    % Get uniformly different values
    idx_z = round( linspace( 1, numel(Zuq), optTree.numSplit ) );
    Zuq_iter = Zuq( idx_z );
    Zuq_iter = unique( Zuq_iter );
else 
    Zuq_iter = Zuq;
end


%% Calculate the splitting criterion
n_split = numel( Zuq_iter );
split_crit = nan(n_split,1);
z_split_id = Z > Zuq_iter(1);
z_split_id_est = Z > Zuq_iter(1);
bw_j = struct;
bw_j.h_l = optTree.bw;
bw_j.h_r = optTree.bw;
bw_j.b_l = optTree.bw_b;
bw_j.b_r = optTree.bw_b;
sq_l = get_np_split_criterion(Y, X, bw_j, c, z_split_id, W, cl, optTree, N_tr, N_est, z_split_id_est );
sq_r = get_np_split_criterion(Y, X, bw_j, c, ~z_split_id, W, cl, optTree, N_tr, N_est, ~z_split_id_est );
split_crit( 1 ) = sq_l + sq_r;
n_tr_i = sum( Ztr <= Zuq_iter( 1 ) & w_tr_log );
n_co_i = sum( Zco <= Zuq_iter( 1 ) & w_co_log );

for i = 2:n_split

    % Check for bucketing
    n_tr_ip = sum( Ztr <= Zuq_iter( i ) & w_tr_log );
    n_co_ip = sum( Zco <= Zuq_iter( i ) & w_co_log );
    n_tr_new = n_tr_ip - n_tr_i;
    n_co_new = n_co_ip - n_co_i;
    if ( n_tr_new >= optTree.obsBucket && n_co_new >= optTree.obsBucket )
        z_split_id = Z > Zuq_iter( i );
        sq_l = get_np_split_criterion(Y, X, bw_j, c, z_split_id, W, cl, optTree, N_tr, N_est, z_split_id_est );
        sq_r = get_np_split_criterion(Y, X, bw_j, c, ~z_split_id, W, cl, optTree, N_tr, N_est, ~z_split_id_est );
        split_crit( i ) = sq_l + sq_r;
        % Update the last valid obs
        n_tr_i = n_tr_ip;
        n_co_i = n_co_ip;
    end

    
    

end

%% Use buckets for more robust criterion selection:
%   Athey-Imbens (2016)
% Check for possible NaN values
logValid = ~isnan( split_crit );

% If there are any candidates
if any( logValid )

    [ ~ , minID ] = min( split_crit );
    
    split = Zuq( minID );

else
    split = NaN;
end

end