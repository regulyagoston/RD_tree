%%
clear all
clc
%
warning('off', 'all');

%% Insert your path to save the results
path = ''

%%
types = {'rdrobust-M3-1tr','rdrobust-M1-2tr','rdrobust-M2-cntstr'};

sigma2_eps_2 = 0.05^2;
MC     = 1000;
N      = [ 1000 , 5000 , 10000, 50000 ];
SE1    = false;

%% Paralell processing is used!


for i = 1 : nT
    for  k = 1 : nS
        disp([i,k])
        mc_tree_fRDD_np_unified( types{ i } , sigma2_eps , N( k ) , SE1, MC, path );
        disp('Unified-tree: done!')
        mc_forest_fRDD_unified( types{ i } , sigma2_eps , N( k ) , SE1 , MC, path );
        disp('Unified forest: done!')
    end
end
