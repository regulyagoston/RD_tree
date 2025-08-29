%% Testing sharp RDD

function mc_forest_fRDD_unified( type , sigmaY , N , SE1 , MC, path )


% Training and Estimation sample
[ Y , Z , X, W ] = fRDDtreeDGPs( type , MC , N , sigmaY , 777 );
N = size( Y , 1 );
% Test Sample
N_te = 10000;
[ ~ , Z_te , ~ , ~, cutoff , chiY ] = fRDDtreeDGPs( type , MC , N_te , sigmaY , 778 );


%% CART options - Simple estimation
optTree                 = optCART;
optTree.maxNodes        = 50;
optTree.maxLevel        = 10;
optTree.maxLeaves       = 15;
optTree.minObs          = 50;
optTree.cp              = 0.0001;%-Inf;
optTree.maxIterGrow     = Inf;
optTree.numKfold        = 10;
optTree.numSplit        = 200;
optTree.honest          = true;
optTree.CV1SE           = false;
optTree.seed            = 1;
optTree.splitType       = 'causal-tree-npRDD-uni';
optTree.type            = 'CLATE';
optTree.model           = 'local';
%optTree.eval_level      = 'leaf-by-leaf';
optTree.criteria        = 'RDD';
optTree.varEstimator    = 'hce-1';%'simple';%
optTree.criteria_weight = 0.5;
optTree.obsBucket       = 4;
optTree.maxBucket       = Inf;
optTree.orderPolinomial = 1;
optTree.kernel          = 'Triang';
optTree.cv_rm_outlier   = true;
optTree.bw_type = 'user-defined';
%optTree.bwselect = 'cerrd';
optTree.forest = true;
optTree.bootstrap_w_replacement = true;
setCritFnc( optTree );

num_tree = 200;
optTree.mtry = get_num_vars_forest( size(Z,2), [] );

% Do not optimize h, take the average from the individual trees
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
filename = horzcat( path, '/fRDD_np_uni_' ,...
    type , '_' , sStr , '_' , nStr , cv_str );
load( filename );
optTree.bw = mean( opt_bw_mc(~isnan(opt_bw_mc)) );
optTree.bw_b = mean( opt_bw_b_mc(~isnan(opt_bw_b_mc)) );


%% MC simulation
numLeaf   = NaN( MC , 1 );
coverage  = NaN( MC , 1 );
bias      = NaN( MC, 1 );
inf_mse   = numLeaf;
%finalTree = cell( MC , 1 );
tau_pred  = cell( MC , 1 );
tau_se_pred  = cell( MC , 1 );


parfor mc = 1 : MC
    
    %% Select sample to use
    obj_sample = samples( optTree , Y( : , mc ) , Z , 'index4TrEst' , [ 1 : N / 2 ; N / 2 + 1 : N ]' ,...
                              'X' , X , 'cutoff' , cutoff, 'W', W );    
    %% Grow Forest
    [ forest , ~, inbag ] = runForest( obj_sample , optTree , num_tree );
    %finalTree{ mc } = forest;
    %% Evaluate infeasible MSE: get the 'predicted' treatment values for given features    
    %% Collect estimates based on true Tree
    [ tau_te_true , tau_var_true ] = predict_forest( forest , inbag , Z_te );
    inf_mse( mc ) = mean( ( chiY - tau_te_true ) .^ 2 );
    tau_pred{ mc }    = tau_te_true;
    tau_se_pred{ mc } = sqrt( tau_var_true ); 
    
    coverage(mc) = get_coverage_rf( tau_pred{ mc } , tau_se_pred{ mc } , chiY , 0.95);
    bias(mc) = mean( chiY - tau_pred{ mc } );
    
end

%% Outputs




filename = horzcat( path, '/fRDD_forest_uni_' , type , '_' , sStr , '_' , nStr , cv_str );
save( filename , 'inf_mse','tau_pred','tau_se_pred','coverage','bias');

end
