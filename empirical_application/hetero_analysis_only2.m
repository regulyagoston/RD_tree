%% Heterogeneous treatment effects in Pop-Elches
clear all
clc

%% Add your Path
path0 = '';
df_all = readtable( horzcat( path0 , 'data/Pop_Elches/dataAER1_R_bw01.csv' ) );
df = df_all( df_all.dzag >= -1.00 & df_all.dzag <= 1.00 & df_all.dzag~=0  , : ); % To speed up the estimation
clear df_all

%% CART options - Simple estimation
optTree                 = optCART;
optTree.maxNodes        = 50;
optTree.maxLevel        = 10;
optTree.maxLeaves       = 8;
optTree.minObs          = 50;
optTree.cp              = 0.0001;
optTree.maxIterGrow     = Inf;
optTree.numKfold        = 10;
optTree.numSplit        = 20;
optTree.honest          = true;
optTree.CV1SE           = false;
optTree.splitType       = 'causal-tree-npRDD-uni';
optTree.type            = 'CATE';
optTree.model           = 'local';
%optTree.eval_level      = 'leaf-by-leaf';
optTree.criteria        = 'RDD';
optTree.varEstimator    = 'hce-1';
optTree.criteria_weight = 0.5;
optTree.obsBucket       = 4;
optTree.maxBucket       = Inf;
optTree.orderPolinomial = 1;
optTree.kernel          = 'Triang';
optTree.cv_rm_outlier   = true;
optTree.bw_type = 'unified';
optTree.bwselect = 'cerrd';
optTree.paralell = true;
setCritFnc( optTree );



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
%Y_bct = df.bct - Dn * b_2fe( 4 : end );
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


%% Unified bandwidth

optTree = copyWithChange(optTree,'bw_type', 'unified' );
num_bw = 8;

[ est_agus , full_agus , optB_a , s_oOS_a, opt_bw_agus, opt_bw_b_agus ] = runTree_uni( samp_agus , optTree, 'grid', num_bw );
[ est_bct  , full_bct  , optB_t , s_oOS_t, opt_bw_bct, opt_bw_b_bct ] = runTree_uni( samp_bct , optTree, 'grid', num_bw );
[ est_bcg  , full_bcg  , optB_g , s_oOS_g, opt_bw_bcg, opt_bw_b_bcg ] = runTree_uni( samp_bcg , optTree, 'grid', num_bw );


%% Get treatment effect along one feature
neval = 25;
% AGUS
f1_a = pop_elches_one_feature( est_agus , samp_agus , 'agus' , 1 , neval );
%%
f2_a = pop_elches_one_feature( est_bct , samp_bct , 'bct' , 1 , neval );
%%
f3_a = pop_elches_one_feature( est_bcg , samp_bcg , 'bcg' , 1 , neval );

%% Nschool
f1_s = pop_elches_one_feature( est_agus , samp_agus , 'agus' , 2 , 25 );
%%
f2_s = pop_elches_one_feature( est_bct , samp_bct , 'bct' , 2 , 25 );
%%
f3_s = pop_elches_one_feature( est_bcg , samp_bcg , 'bcg' , 2 , 25 );