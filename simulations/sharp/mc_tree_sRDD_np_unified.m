%% Testing sharp RDD

function mc_tree_sRDD_np_unified( type , sigmaY , N , SE1, MC, path )

% Training and Estimation sample
[ Y , Z , X ] = sRDDtreeDGPs( type , MC , N , sigmaY , 777 );
N = size( Y , 1 );
% Test Sample
N_te = 10000;
[ ~ , Z_te , ~ , cutoff , chiY , ~ , predVals , trueTree ] = sRDDtreeDGPs( type , MC , N_te , sigmaY , 778 );


%% CART options - Simple estimation
optTree                 = optCART;
optTree.maxNodes        = 50;
optTree.maxLevel        = 10;
optTree.maxLeaves       = 15;
optTree.minObs          = 50;
optTree.cp              = 0.0001;
optTree.maxIterGrow     = Inf;
optTree.numKfold        = 10;
optTree.numSplit        = 200;
optTree.honest          = true;
optTree.CV1SE           = false;
optTree.seed            = 1;
optTree.splitType       = 'causal-tree-npRDD-uni';
optTree.type            = 'CATE';
optTree.model           = 'local';
optTree.criteria        = 'RDD';
optTree.varEstimator    = 'hce-1';%'simple';%
optTree.criteria_weight = 0.5;
optTree.obsBucket       = 4;
optTree.maxBucket       = Inf;
optTree.orderPolinomial = 1;
optTree.kernel          = 'Triang';
optTree.cv_rm_outlier   = true;
optTree.bw_type = 'unified';
optTree.bwselect = 'cerrd';
setCritFnc( optTree );


%% MC simulation
numLeaf   = NaN( MC , 1 );
coef_se   = NaN( MC , 1 );
dgp_found = NaN( MC , 1 );
coverage  = NaN( MC , 1 );
bias  = NaN( MC , 1 );
inf_mse   = numLeaf;
tau_pred  = cell( MC , 1 );
tau_se_pred  = cell( MC , 1 );
opt_bw_mc = numLeaf;
opt_bw_b_mc = numLeaf;
opt_beta_mc = numLeaf;


parfor mc = 1 : MC
    optTree_i = copyWithChange(optTree,'num_bw', 10);
    %% Select sample to use
    obj_sample = samples( optTree_i , Y( : , mc ) , Z , 'index4TrEst' , [ 1 : N / 2 ; N / 2 + 1 : N ]' ,...
                              'X' , X , 'cutoff' , cutoff );    
    %% Grow and find optimal tree
    % Cross-validation
    [ opt_beta , opt_bw, opt_bw_b, s_oOS, bw, beta ] = cross_validate_bws_grid( obj_sample , optTree_i, [] );
    optTree_i = copyWithChange(optTree_i,'bw', opt_bw);
    optTree_i = copyWithChange(optTree_i,'bw_b', opt_bw_b);
    opt_bw_mc( mc ) = opt_bw;
    opt_beta_mc( mc ) = opt_beta;
    opt_bw_b_mc( mc ) = opt_bw_b;
    obj_tree = growRDDtree( obj_sample , optTree_i );
    finalTree_tr = pruning( obj_tree , obj_sample , optTree_i , opt_beta );

    % Number of leaves
    numLeaf( mc ) = numel( findleaves( finalTree_tr ) );
    
    %% Estimate proper CATE
    finalTree = estimate_tree_est( finalTree_tr , obj_sample , optTree_i );
    
    %% Evaluate infeasible MSE: get the 'predicted' treatment values for given features
    tau_te = predict( finalTree, Z_te , 'tau_j' );
    % Infeasible MSE
    inf_mse( mc ) = mean( ( chiY - tau_te ) .^ 2 );
    
    %% Check whether the true tree is found
    dgp_found( mc ) = issameTree( finalTree , trueTree , 0.5 );
    %% Collect estimates based on true Tree
    tau_pred{ mc }    = predict( finalTree , predVals , 'tau_j' );
    tau_se_pred{ mc } = predict( finalTree , predVals , 'tau_se_j'); 
    
    %% Try Athey-Imbens variance estimator
    coef_var = get_var_tree( finalTree , obj_sample , false );
    coef_se( mc ) = sqrt( coef_var );
    
    %% Coverage for the test sample
    [ coverage( mc ), bias( mc ) ] = get_coverage( finalTree, Z_te , chiY , 0.95 );
    
end

%% Outputs

if sigmaY < 1
    aux_s2 = num2str( sigmaY );
    sStr = horzcat('0',aux_s2( 3 : end ) );
else
    sStr = num2str( sigmaY );
end
if SE1
    cv_str = '_1SE';
else
    cv_str = '_noSE';
end

nStr = num2str( N );
filename = horzcat( path, '/sRDD_np_uni_' , type , '_' , sStr , '_' , nStr , cv_str );
save( filename , 'inf_mse','numLeaf','dgp_found', 'tau_pred','tau_se_pred','coef_se','coverage','opt_bw_mc','bias','opt_bw_b_mc');

end
