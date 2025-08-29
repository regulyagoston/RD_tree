%% Get Monte-Carlo simulation variables for sharp `tree' from Reguly (2021)
%
% Input:
%   type:
%   Calonico et al.:
%       'rdrobust-M1-2tr','rdrobust-M2-cntstr','rdrobust-M3-1tr'
%   M: number of MC iteration
%   N: number of obs
%   sigma2_eps: variance of disturbance term (N(0,sigma^2_eps))
%   seed: for random number generation
%
% Output:
%   Y - outcome variable (Nx1)
%   Z - features (one or two variable)
%   W - treatment status (Nx1)
%   chi - true conditional (average) treatment function
%   eta - true conditional mean function for Y apart getting treatment or not
%   predValues - true predicted tau values for evaluating estimated values
%   trueTree - true tree structure if exists


function [ Y , Z , X , c , chi , eta , predValues , trueTree ] = sRDDtreeDGPs( type , M , N , sigma2_eps , seed )

rng( seed );

% Definition of getting treatment
X = rand( [ N , 1 ] ) * 2 - 1;
c = 0;
W = X >= c;
predValues = NaN;
nPred = 100;
sigmaY = sqrt( sigma2_eps );

% Initialize true tree
trueTree = tree;

switch type   
%% Calonico et al simulations        
    case {'rdrobust-M1-2tr','rdrobust-M2-cntstr','rdrobust-M3-1tr'}
        % Definition of getting treatment
        X = random( 'beta' , 2 , 4 , [ 2*N , 1 ] ) * 2 - 1;
        dropX = X < -0.99 | X > 0.99;
        X = X( ~dropX );
        X = X( 1 : N );
        c = 0;
        W = X >= c;
        eps = randn( [ N , M ] ) .* sigmaY;
        rootNode = nodeProp( [] , [] , [] , 0 , true( 1 , 50 ) );
        trueTree = addnode( trueTree , 0 , rootNode );
        if any( strcmp( type , {'rdrobust-M1-2tr','rdrobust-M3-1tr'} ) )
            Z_pol    = rand( N , 1 ) > 0.5;
            n_states = floor( N ./ 50 );
            Z_states = zeros( N , 49 );
            for i = 1 : 48
                Z_states( ( i - 1 ) * n_states + 1 : i * n_states , i ) = 1;
            end
            Z_states( i * n_states + 1 : end ,  end ) = 1;
            Z = [ Z_pol , Z_states( randperm( N ) , : ) ];
        end
        if strcmp( type , 'rdrobust-M1-2tr' )
            % Democrats and republicans has different incumbancy advantage
            chi  = ( 0.56 - 0.48 ) .* Z( : , 1 ) + ( 0.50 - 0.48 ) .* ( 1 - Z( : , 1 ) );
            eta1_1 = 0.48 + 1.27 .* X + 7.18 .* X.^2 + 20.21 .* X.^3 + 21.54 .* X.^4 + 7.33 .* X.^5;
            eta1_2 = 0.48 + 2.35 .* X + 8.18 .* X.^2 + 22.21 .* X.^3 + 24.14 .* X.^4 + 8.33 .* X.^5;
            eta2_1 = 0.48 + 0.84 .* X - 3.00 .* X.^2 +  7.99 .* X.^3 -  9.01 .* X.^4 + 3.56 .* X.^5;
            eta2_2 = 0.48 + 1.21 .* X - 2.90 .* X.^2 +  6.99 .* X.^3 - 10.01 .* X.^4 + 4.56 .* X.^5;
            Y = repmat( ( eta1_1 .* Z( : , 1 ) + eta1_2 .* ( 1 - Z( : , 1 ) ) ) .* ( 1 - W ) +...
                        ( eta2_1 .* Z( : , 1 ) + eta2_2 .* ( 1 - Z( : , 1 ) ) ) .* W + W .* chi , [ 1 , M ] ) + eps;
            if nPred*2 > N
                nPred = floor( N ./ 2 );
            end
            eta = [ eta1_1 , eta1_2 , eta2_1 , eta2_2 ];
            ids = createPartition( nPred * 2 , 1 , seed );
            predValues = [ [ zeros( nPred , 1 ) ; ones( nPred , 1 ) ] , Z( ids , 2 : end ) ];
            % Add the two candidate values
            node_1 = nodeProp( 1 , 0 , true  , 0 , true( 1 , 50 ) );
            node_2 = nodeProp( 1 , 0 , false , 0 , true( 1 , 50 ) );
            trueTree = addnode( trueTree , 1 , node_1 );
            trueTree = addnode( trueTree , 1 , node_2 );
        elseif strcmp( type , 'rdrobust-M2-cntstr' )
            % Different Continents
            n_cont = floor( N ./  6 );
            Z_cont = zeros( N , 5 );
            for i = 1 : 4
                Z_cont( ( i - 1 ) * n_cont + 1 : i * n_cont , i ) = 1;
            end
            Z_cont( i * n_cont + 1 : end , end ) = 1;
            Z_age = rand( N , 1  ) * 4;%  + 5
            Z = [ Z_age , Z_cont( randperm( N ) , : ) ];
            %chi  = ( 0.26 - 3.71 )
            chi  = -0.45 + 0.5*Z_age - 0.25*Z_age.^2 + 0.1*Z_age.^3; %- 5*Z_age;%
            eta1 = 3.71 +  2.30 .* X +  3.28 .* X.^2 +  1.45 .* X.^3 +  0.23 .* X.^4 + 0.03 .* X.^5;
            eta2 = 3.71 + 18.49 .* X - 54.81 .* X.^2 + 74.30 .* X.^3 - 45.02 .* X.^4 + 9.83 .* X.^5;
            Y = repmat( eta1 .* ( 1 - W ) + eta2 .* W + W .* chi , [ 1 , M ] ) + eps;
            if nPred*5 > N
                nPred = floor( N ./ 5 );
            end
            eta = [ eta1 , eta2 ];
            ids = createPartition( nPred , 1 , seed );
            predValues = [ sort( repmat( [ 5 ; 6 ; 7 ; 8 ; 9 ] , [ nPred ./ 5 , 1 ] ) ) , ...
                            Z( ids , 2 : end ) ];
            % add rootenode
            rootNode = nodeProp( [] , [] , [] , 0 , true( 1 , 6 ) );
            trueTree = addnode( trueTree , 0 , rootNode );
        elseif strcmp( type , 'rdrobust-M3-1tr' )
            chi  = ( 0.52 - 0.48 ) + zeros( N , 1 );
            eta1 = 0.48 + 1.27 .* X - 0.5 .* 7.18 .* X.^2 + 0.7 .* 20.21 .* X.^3 + 1.1 .* 21.54 .* X.^4 + 1.5 .* 7.33 .* X.^5;
            eta2 = 0.48 + 0.84 .* X - 0.1 .* 3.00 .* X.^2 - 0.3 .*  7.99 .* X.^3 - 0.1 .*  9.01 .* X.^4 + 3.56 .* X.^5;
            Y = repmat( eta1 .* ( 1 - W ) + eta2 .* W + W .* chi , [ 1 , M ] ) + eps;
            eta = [ eta1 , eta2 ];
            ids = createPartition( nPred , 1 , seed );
            predValues = Z( ids , : );
        end
    otherwise
        error('No such simulation setup')
end

end