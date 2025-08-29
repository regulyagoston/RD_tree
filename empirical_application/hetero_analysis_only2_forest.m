%% Heterogeneous treatment effects in Pop-Elches
clear all
clc

%% Path
path0 = '';
df_all = readtable( horzcat( path0 , 'data/Pop_Elches/dataAER1_R_bw01.csv' ) );
df = df_all( df_all.dzag >= -1.00 & df_all.dzag <= 1.00 & df_all.dzag~=0  , : ); % 
clear df_all

%% CART options - Simple estimation
optTree                 = optCART;
optTree.maxNodes        = 50;
optTree.maxLevel        = 10;
optTree.maxLeaves       = 15;
optTree.minObs          = 50;
optTree.cp              = 0.0001;%-Inf;
optTree.maxIterGrow     = Inf;
optTree.numKfold        = 10;
optTree.numSplit        = 20;
optTree.honest          = true;
optTree.CV1SE           = false;
optTree.seed            = 1;
optTree.splitType       = 'causal-tree-npRDD-uni';
optTree.type            = 'CATE';
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


%% Creating feature variables
% Number of schools in town
[ d_nusua , k_nusua ] = create_dummies( df.nusua );
% Creating the features - only variables as in Pop-Elches Table 5
Z = [ df.agus , df.nusua , d_nusua ];

%% Set up the samples

% Running variable
X = df.dzag;
% Clusters
cl = df.sid2;

% School lvlv avg transition score
Y_agus = df.agus_dm;

% BA prob
Y_bct = df.bct_dm;

% BA grade
lN = ~isnan( df.bcg_dm );
Y_bcg = df.bcg_dm( lN );

% Training and estimation samples
N_12 = numel( df.agus );
idxTrEst_12 = createPartition( N_12 , 2 , 777 );
N_3 = sum( lN );
idxTrEst_3 = createPartition( N_3 , 2 , 777 );

% Samples
samp_agus = samples( optTree , Y_agus , Z , 'index4TrEst' , idxTrEst_12 ,...
                     'X' , X , 'cutoff' , 0 , 'cluster' , cl );
samp_bct = samples( optTree , Y_bct , Z , 'index4TrEst' , idxTrEst_12 ,...
                     'X' , X , 'cutoff' , 0 , 'cluster' , cl );
samp_bcg = samples( optTree , Y_bcg , Z( lN , : ) , 'index4TrEst' , idxTrEst_3 ,...
                     'X' , X( lN ) , 'cutoff' , 0 , 'cluster' , cl( lN ) );

%% Re-Run with unified bandwidth

% Values from individual trees
optTree = copyWithChange(optTree,'bw', 0.1481 ); % bw value comes from optimal tree estimate
optTree = copyWithChange(optTree,'bw_b', 0.1481 );
[ est_agus , ~ , inbag_agus ] = runForest( samp_agus , optTree, num_tree );
save('pop_elches_forest_agus.mat','est_agus','-v7');
save('pop_elches_forest_agus_inbag.mat','inbag_agus');
save('pop_elches_forest_agus.mat','est_agus');
optTree = copyWithChange(optTree,'bw', 0.0836 ); % bw value comes from optimal tree estimate
optTree = copyWithChange(optTree,'bw_b', 0.0836 );
[ est_bct  , ~ , inbag_bct ] = runForest( samp_agus , optTree, num_tree );
save('pop_elches_forest_bct.mat','est_bct','inbag_bct');
optTree = copyWithChange(optTree,'bw', 0.2247 ); % bw value comes from optimal tree estimate
optTree = copyWithChange(optTree,'bw_b', 0.2247 );
[ est_bcg  , ~ , inbag_bcg ] = runForest( samp_agus , optTree, num_tree );


%% Marginal treatments
[ tau_agus , tau_var_agus ] = predict_forest( est_agus , inbag_agus , samp_agus.Z_est );

% I have added maually these estimates to the tree estimate using the code for: hetero_analysis_only2.m and the function `pop_elches_one_feature`