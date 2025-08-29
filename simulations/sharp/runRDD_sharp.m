%%
clear all
clc
%
% To do: rerun rdrobust-M2-cntstr for all, rdrobust-M1-2tr for N=10,000
% with NP and from 'cnts-hetero-cntstr' the remaining with N=50,000
warning('off', 'all');

%% DEFINE YOUR PATH TO SAVE RESULTS
path = '';


types = {'rdrobust-M3-1tr','rdrobust-M1-2tr','rdrobust-M2-cntstr'};

sigma2_eps_2 = 0.05^2;
MC     = 1000;
N      = [1000,5000,10000,50000];
SE1    = false;
mc=1;


for i = 1 : nT
    for  k = 1 : nS
        disp([i,k])
        mc_tree_sRDD_np_unified( types{ i } , sigma2_eps , N( k ) , SE1, MC, path );
        disp('Unified tree: done!')
        mc_forest_sRDD_unified( types{ i } , sigma2_eps , N( k ) , SE1 , MC, path );
        disp('Unified forest: done!')
    end
end

