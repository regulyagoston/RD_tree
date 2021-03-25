%% Randomly partition 1 , ... , N sample into k subsamples:
% Agoston Reguly (2018)
%
% Input:
%   N - 1x1 positive integer for sample size
%   k - 1x1 positive integer for the number of subsamples
%
% Output:
%   kIdx - kmaxN x k index matrix, containing the indexes for each
%           subsamples
%   kmaxN - 1x1 positive integer, maximal sample size for the subsamples
%
% Last modified: 19.06.2018

function [ kIdx , kmaxN ] = createPartition( N , k , seed )

if numel( N ) ~= 1 && N < 1 && mod( N , 1 ) ~= 0
    error('N must be a positive integer for partitioning')
end
if numel( k ) ~= 1 && k < 1 && mod( k , 1 ) ~= 0
   error('k must be a positive integer!') 
end

% Number of maximal elements in each fold
kmaxN = ceil( N / k );
% Number of maximal elements in each fold
kminN = floor( N / k );
% Create the new index matrix
kIdx = NaN( kmaxN , k );
% Randomly choose the indexes
rng( seed );
rIdx = randperm( N , N );
if kmaxN == kminN
    kIdx = reshape( rIdx , kmaxN , k );
else
    % Reshape indexes to all slots, excluding the last indexes
    kIdx( 1 : end - 1 , : ) = reshape( rIdx( 1 : end - ( N - k * kminN  ) ) , kminN , k );
    % Fill the remaining indexes
    for i = 1 : k * kmaxN - N
        kIdx( end , i ) = rIdx( end - ( k * kmaxN - N ) + i );
    end
    % Double the last indexes, where there is less observation
    kIdx( end , i + 1 : end ) = kIdx( end - 1 , i + 1 : end );
end

if any( all( isnan( kIdx ) , 2 ) )
    kIdx = kIdx( ~all( isnan( kIdx ) , 2 ) , : );
end

end