%% Example for Regression Discontinuity Tree
%   Agoston Reguly (2025)
clear all
clc

%% Generating simple RD following Calanico et al (2014)
rng(777)
N  = 10000;
sigmaY = 0.05^2;
% Definition of getting treatment
X = random( 'beta' , 2 , 4 , [ 2*N , 1 ] ) * 2 - 1;
dropX = X < -0.99 | X > 0.99;
X = X( ~dropX );
X = X( 1 : N );
c = 0;
W = X >= c;
eps = randn( [ N , 1 ] ) .* sigmaY;
Z_pol    = rand( N , 1 ) > 0.5;
n_states = floor( N ./ 50 );
Z_states = zeros( N , 49 );
for i = 1 : 48
    Z_states( ( i - 1 ) * n_states + 1 : i * n_states , i ) = 1;
end
Z_states( i * n_states + 1 : end ,  end ) = 1;
Z = [ Z_pol , Z_states( randperm( N ) , : ) ];
chi  = ( 0.56 - 0.48 ) .* Z( : , 1 ) + ( 0.50 - 0.48 ) .* ( 1 - Z( : , 1 ) );
eta1_1 = 0.48 + 1.27 .* X + 7.18 .* X.^2 + 20.21 .* X.^3 + 21.54 .* X.^4 + 7.33 .* X.^5;
eta1_2 = 0.48 + 2.35 .* X + 8.18 .* X.^2 + 22.21 .* X.^3 + 24.14 .* X.^4 + 8.33 .* X.^5;
eta2_1 = 0.48 + 0.84 .* X - 3.00 .* X.^2 +  7.99 .* X.^3 -  9.01 .* X.^4 + 3.56 .* X.^5;
eta2_2 = 0.48 + 1.21 .* X - 2.90 .* X.^2 +  6.99 .* X.^3 - 10.01 .* X.^4 + 4.56 .* X.^5;
Y = ( eta1_1 .* Z( : , 1 ) + eta1_2 .* ( 1 - Z( : , 1 ) ) ) .* ( 1 - W ) +...
    ( eta2_1 .* Z( : , 1 ) + eta2_2 .* ( 1 - Z( : , 1 ) ) ) .* W + W .* chi + eps;
% Scatter plot to see the design
% scatter( X , Y )


%% CART options - Simple estimation%% CART options - Simple estimation
optTree                 = optCART;
optTree.maxNodes        = 50; % Max number of nodes in the tree
optTree.maxLevel        = 10; % Maximum depth of the tree
optTree.maxLeaves       = 15; % Maximum number of leaves
optTree.minObs          = 50; % Minimum number of effective observations within each leaf
optTree.cp              = 0.0001; % Criterion improvement requirement for leaf to be valid
optTree.maxIterGrow     = Inf; % Maximum number of iteration by the algorithm
optTree.numKfold        = 10;  % Number of folds during cross-validation
optTree.numSplit        = Inf; % Maximum number of splits for each feature
optTree.honest          = true; % Honest sample splitting
optTree.CV1SE           = false; % For cross-validation use 1SE rule or not
optTree.seed            = 1; % Seed
optTree.splitType       = 'causal-tree-npRDD-uni'; % Split criterion -- use unified bandwidth with nonparametric estimator
optTree.type            = 'CATE'; % CATE or CLATE
optTree.model           = 'local'; % local polynomial model
optTree.criteria        = 'RDD';   % Criterion is for RDD
optTree.varEstimator    = 'hce-1'; % Type of variance estimator
optTree.obsBucket       = 4; % Minimum number of effective observations required for valid splitting value
optTree.orderPolinomial = 1; % Order of polynomial
optTree.kernel          = 'Triang'; % Kernel type
optTree.cv_rm_outlier   = true; % Exclude outliers during cross-validation or not
optTree.bw_type = 'unified'; % Type of bandwidth
optTree.bwselect = 'cerrd';  % Initial bandwidth is given by CE or MSE optimal bw
setCritFnc( optTree );

%% Set-up the sample object
% Indexes for the training sample and for the estimation sample:
%   here just simply split the sample in the middle, but can use any
%   sampling tecnique.
if mod( N / 2 , 2 ) == 0
    indTrEst = [ 1 : N / 2 ; N / 2 + 1 : N ]';
else
    indTrEst = [ 1 : ( N + 1 ) / 2 ; [ ( N + 1 ) / 2 + 1 : N , N ] ]';
end
% Set up the sample
obj_sample = samples( optTree , Y , Z , 'index4TrEst' , indTrEst , 'X' , X , 'cutoff' , c );
% can use 'cluster' and the variable to estimate clustered SE

%% Find and estimate the optimal tree
bw_search = 'grid';
num_bw = 20;
[ final_tree , large_tree, opt_gamma, opt_bw, emse_crit, cand_bw ] = runTree_uni( obj_sample , optTree, bw_search, num_bw );


% Outputs:
%  final_tree: estimated optimal RD tree
% not necessary outputs:
%  large_tree: large tree estimated on training sample
%  opt_gamma: optimal pruning parameter
%  opt_bw: optimal bandwidth parameter
%  emse_crit: averaged cross-validation criteria for each bandwidth
%  parameters (with optimal prining parameter)
%  cand_bw: candidate bandwidth parameters

%%
% Visualize the tree:
tostring( final_tree )

% Print optimal parameters
sprintf(horzcat('Optimal bandwidth parameter:', num2str( opt_bw ) ) )
sprintf(horzcat('Optimal pruning parameter:', num2str( opt_gamma ) ) )