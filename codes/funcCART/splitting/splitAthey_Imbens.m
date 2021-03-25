%% Split function Athey-Imbens - using minimum number of buckets
% 
% Input:
%   Z_i - Nx1 vector
%       feature to be splitted
%   numSplit - scalar
%       number of maximum splits
%   W_i - Nx1 vector
%       treatment indicator (W_i and Z_i must be matched)
%   b0 - scalar
%       target number of observations in a bucket
%   maxBucket - scalar
%       maximum nuber of buckets
%   minSize - scalar
%       minimum number of buckets
%   
% Example:
%   Z_i = rand( 100 , 1 );
%   numSplit = 20;
%   Wi = rand( 100 , 1 ) > 0.5;
%   b0 = 4;
%   maxBucket = 50;
%   minSize   = 10;


function split = splitAthey_Imbens( Z_i , numSplit , Wi , b0  , maxBucket , minSize )

Z_t = sort( Z_i( Wi == 1 ) );
nZt = numel( Z_t );
Z_c = sort( Z_i( Wi == 0 ) );
nZc = numel( Z_c );

% Get the unique values which are potential splits
Z_u = unique( Z_i );
nU  = numel( Z_u );

% Count each unique observations
num_obs_treat = NaN( 1 , nU );
num_obs_contr = NaN( 1 , nU );
for i = 1 : nU
    num_obs_treat( i ) = sum( Z_u( i ) == Z_t );
    num_obs_contr( i ) = sum( Z_u( i ) == Z_c );
end

% Cummulative observations
co_treat = cumsum( num_obs_treat );
co_contr = cumsum( num_obs_contr );

% Set the lower and the upper split values to have enough values
aux_1 = Z_u( co_treat >= minSize );
aux_2 = Z_u( co_contr >= minSize );
aux_3 = Z_u( co_treat >= nZt - minSize );
aux_4 = Z_u( co_contr >= nZc - minSize );
split = [ max( aux_1( 1 ) , aux_2( 1 ) ) , min( aux_3( 1 ) , aux_4( 1 ) ) ];
% Return if there is not enough observations return
if numel( split ) < 2
    return;
end
% If there is an overlap, choose the first split value and return
if split( 2 ) < split( 1 )
    split = split( 1 );
    return;
end
% Adjust the number of splits
numSplit = numSplit - 2;
% Remove values
log_candidates   = Z_u > split( 1 ) & Z_u < split( 2 );
split_candidates = Z_u( log_candidates );
% If not enough observations
nSc = numel( split_candidates );
if nSc <= numSplit
    split = [ split , split_candidates' ];
    return;
end
% Set the target number of observations in a bucket accordingly
if ceil( nSc ./ b0 ) < numSplit
    b0 = max( floor( nSc ./ numSplit ) , 1 );
end
if ceil( nSc ./ b0 ) > maxBucket
    b0 = max( floor( nSc ./ maxBucket ) , 1 );
end
    
% Adjust the cummulative sum
aux_11 = find( Z_u <= split( 1 ) , 1 , 'last' );
co_treat_2 = co_treat( log_candidates ) - co_treat( aux_11 );
co_contr_2 = co_contr( log_candidates ) - co_contr( aux_11 );

%% Find the optimal split points
maxIter = 100;
b       = NaN( maxIter , 1 );
b( 1 )  = b0;
max_obs = min( co_treat_2( end ) , co_contr_2( end ) );
for i = 1 : maxIter
    split_c = [];
    co_treat_3 = co_treat_2;
    co_contr_3 = co_contr_2;
    for j = 1 : ceil( max_obs ./ b( i ) )
        logSplit = co_treat_3 <= b( i ) | co_contr_3 <= b( i );
        if ~all( logSplit )
            % Get the last valid split candidate value
            aux_1    = split_candidates( logSplit );
            split_c  = [ split_c , aux_1( end ) ];
            % Set those observations to zeros and substract the last value
            aux_2    = co_treat_3( logSplit );
            co_treat_3(  logSplit ) = 0;
            co_treat_3( ~logSplit ) = co_treat_3( ~logSplit ) - aux_2( end );
            % Set those observations to zeros and substract the last value
            aux_3    = co_contr_3( logSplit );
            co_contr_3(  logSplit ) = 0;
            co_contr_3( ~logSplit ) = co_contr_3( ~logSplit ) - aux_3( end );
        else
            break;
        end
    end
    % We use the minimum of these two possible buckets
    nB = numel( split_c );
    % Case of too few splits -> decrease the number of observations in the buckets
    if nB < numSplit
        b( i + 1 ) = b( i ) - 1;%ceil( ( numSplit - nB ) / b( i ) );
        tr_S = false;
    else
        tr_S = true;
    end
    if nB > maxBucket
        b( i + 1 ) = b( i ) + 1;%ceil( ( nB - maxBucket ) / b( i ) );
        tr_B = false;
    else
        tr_B = true;
    end
    if b( i + 1 ) < 1
        b( i + 1 ) = 1;
    end
    if ( tr_S && tr_B ) || b( i ) == 1 || i == maxIter
        split = sort( [ split , split_c ] );
        break;
    end
end

end
