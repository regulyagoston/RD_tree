%% Create CATE plot along one feature


function [ xeval , tau_bs , tau_se_bs ] = cate_plot_onefeature( tree , sample , Zid , neval )

% Get the treatments and the features
[ tau_all , tau_se_all , Z_all ] = get_treatment_n_features( tree , sample );
% Sort along the selected feature
[ Zs , idZ ] = sort( Z_all( : , Zid ) );
tau_sorted = tau_all( idZ );
tau_se_sorted = tau_se_all( idZ );
% Create bin-scatter values (using equal spacing and means
if Zid == 2
    xeval = min( sample.Z_est( : , Zid ) ) : max( sample.Z_est( : , Zid ) );
else
    xeval = linspace( min( sample.Z_est( : , Zid ) ) , max( sample.Z_est( : , Zid ) ) , neval );
end
tau_bs = binscatter( tau_sorted , Zs , xeval );
tau_se_bs = binscatter( tau_se_sorted , Zs , xeval );

% % Create graphs
% f1 = figure('Position', get(0, 'Screensize'));
% if Zid == 1
%     plot( xeval , tau_bs , 'b','linewidth' , 2 );
%     hold on
%     plot( xeval , tau_bs + 2*tau_se_bs , '--', 'color','r','linewidth' , 2 );
%     plot( xeval , tau_bs - 2*tau_se_bs , '--', 'color','r','linewidth' , 2 );
% else
%     errorbar( xeval , tau_bs , 2*tau_se_bs, 'linewidth' , 2 );
% end
% ylabel('Treatment effect')
% box('on');
% grid('on');
% legend('boxoff')
% set(gca,'fontsize',22)

end