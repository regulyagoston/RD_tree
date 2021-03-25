%% Splitting criterion for rpart - prediction, MSE, using means and leaf-by-leaf criterion

function split = split_rpart_classical( Y_leaf , Zi , optTree )


% All the variables to be splitted
Y_temp = Y_leaf;
% Variable is assumed to be continuous -> tricky way to calculate MSE
% 1st sort according to Zi and get observations to satisfy minimum
% number of obs in a leaf
[ Zs , sortID ] = sort( Zi );
Y_temp = Y_temp( sortID , : );
% Demean the outcome
nY     = numel( Y_temp );
Y_sc   = Y_temp - mean( Y_temp );
% Indicating the means by each split
temp    = cumsum( Y_sc( 1 : end ) );
left_nY  = ( 1 : nY )';
right_nY = nY - left_nY;
l_mean  =  temp ./ left_nY;
r_mean  = -temp ./ right_nY;
spQ_i   = -( left_nY .* l_mean .^ 2 + right_nY .* r_mean .^2 ) ./ sum( Y_sc .^ 2 );
% Find valid minimum
[ ~ , mID ] = min( spQ_i( optTree.minObs : end - optTree.minObs ) );
% Get the splitting value
split = Zs( mID + optTree.minObs );

end