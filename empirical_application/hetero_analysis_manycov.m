%% Heterogeneous treatment effects in Pop-Elches
clear all
clc

%% Path
path0 = '';
df_all = readtable( horzcat( path0 , 'data/Pop_Elches/data-AER-7.csv' ) );
df = df_all(df_all.dzag~=0,:);%( df_all.dzag >= -1.00 & df_all.dzag <= 1.00 & df_all.dzag~=0 , : );


%% CART options - Simple estimation
optTree                 = optCART;
optTree.maxNodes        = 50;
optTree.maxLevel        = 10;
optTree.maxLeaves       = 15;
optTree.minObs          = 25;
optTree.cp              = 0.0001;%-Inf;
optTree.maxIterGrow     = Inf;
optTree.numKfold        = 10;
optTree.numSplit        = 400;
optTree.honest          = true;
optTree.CV1SE           = true;
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
optTree.bw_type = 'unified';
optTree.bwselect = 'cerrd';
optTree.paralell = true;
setCritFnc( optTree );

%% Creating feature variables
Znames = {'school tr. score', 'nuum schools' , '2 schools', '3 schools','4 schools',...
          'Female head of household','Age of head of household','Romanian','Hungarian','Gypsy','Other Ethnicity',...
          'Primary educ' ,'Secondary educ' , 'Tertiary educ' , 'Sex of student'  , 'Age of student',...
          'Parent volunteered' , 'Parent paid tutoring','Parent helps HW','Child does HW every day - Parent',...
          'Negative interactions w peers','Child does HW every day - Child', 'HW percieved easy',...
          'Highest certification proportion - teacher','Proportion of novices',...
          'Car','Internet','Phone','computer',...
          'Mom_educ_primary','Mom_educ_secondary','Mom_educ_tertiary'};
[ d_nusua , k_nusua ] = create_dummies( df.nusua );
[ d_tutor , k_tut ] = create_dummies( df.p_tutoring );
ch_rank_peers = str2double( df.ch_rank_peers );%'Rank among peers (1-7)',
exper_rom = str2double( df.experience_Romanian );%'Years of experience'
Z = [ df.agus , df.nusua , d_nusua , ...
      df.head_sex , df.head_age , df.head_nat_romanian , df.head_nat_hungarian , df.head_nat_gypsy , ...
      df.head_nat_other , df.head_educ_primary , df.head_educ_sec , df.head_educ_tertiary , ...
      df.ch_sex , df.ch_age , df.p_d_parent_volunteer , d_tutor( : , 3 ) , df.p_d_homework_help , df.p_d_homework , ...
      df.ch_peers_index_bad , df.ch_d_homework , df.ch_rank_homework_index , ...
      df.didactic_Romanian , df.novice_Romanian , ...
      df.con_car , df.con_internet , df.con_phone , df.con_computer,...
      df.mom_educ_primary,df.mom_educ_secondary,df.mom_educ_tertiary];


lValid = all(~isnan(Z),2); %& l23_schools;
Nv =  sum( lValid );

%% Set up the samples
% Fixed effects
Dn = create_dummies( df.uazY );
% Running variable
X = df.dzag;
% Clusters
cl = df.sid2;

% School lvlv avg transition score
b_1fe = lscov( [ df.dga , df.dzag , df.dzag_after , Dn ] , df.agus );
Y_agus = df.agus - Dn * b_1fe( 4 : end );

% Training and estimation samples
idxTrEst_12 = createPartition( Nv , 2 , 777 );


% Samples
samp_agus = samples( optTree , Y_agus( lValid ) , Z( lValid , : ) , 'index4TrEst' , idxTrEst_12 ,...
                     'X' , X( lValid ) , 'cutoff' , 0 , 'cluster' , cl( lValid ) );

%% Run discontinuity tree
optTree = copyWithChange(optTree,'bw_type', 'unified' );
num_bw = 20;

[ est_agus , full_agus , optB_a , oos, opt_bw, opt_bw_b ] = runTree_uni( samp_agus , optTree, 'grid', num_bw );

%% Create plots for CATE along agus and nusua
optTree.bw = opt_bw;
optTree.bw_b = opt_bw_b;
b1 = optB_a;
finalTree_tr = pruning( full_agus , samp_agus , optTree , b1(end) );
% Do the estimation on the estimation sample
est_tree2 = estimate_tree_est( finalTree_tr , samp_agus , optTree );
tostring(est_tree2)
%%
%tostring( est_agus )
var_id = get_relevant_vars( est_tree2 );
Znames( var_id ) 
a=Znames( var_id );