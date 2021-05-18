%% Randomly partition a sample into k subsamples with N observations
% Agoston Reguly (2018)
%
% Input:
%   X - sample: 1st dimension is the number of observations, 2nd is the
%       number of variables.
%   N - 1x1 number of observations in the sample
%   k - 1x1 positive integer for the number of bootstrap subsamples
%   seed - 1x1 positive integer for set radnom number generator
%
% Output:
%   kIdx - N x k index matrix, containing the indexes for each
%           bootstrap samples
%
% Last modified: 30/03/2021

function kIdx = createBootstrapSampleIdx( X , N , k , seed )

if ~isnumeric( X )
   error( 'X must be a numeric vector/matrix.' ) 
end

if numel( N ) ~= 1 && N < 1 && mod( N , 1 ) ~= 0
    error('N must be a positive integer for number of observations in bootstrap sample!')
end
if numel( k ) ~= 1 && k < 1 && mod( k , 1 ) ~= 0
   error('k must be a positive integer!') 
end

% Number of observations in the matrix
nX = size( X , 1 );

% Randomly choose the indexes
rng( seed );
kIdx = randi( nX , [ N , k ] );


end