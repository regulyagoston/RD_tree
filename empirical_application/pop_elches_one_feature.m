%% Pop-Elches  Replicate Figure A4


function f1 = pop_elches_one_feature( tree , sample , type , Zid , neval )

[ xv_t , tau , tau_se ] = cate_plot_onefeature( tree , sample , Zid , neval );

Zvar = sample.Z_est( : , Zid );
if Zid == 1
    xval = linspace( min( Zvar ) , max( Zvar ) , neval );
else
    xval = unique(Zvar);
    neval = numel( xval );
end
if Zid == 1
    pop_vals = NaN( neval , 3 );
    pop_vals_se = pop_vals;
    if strcmp( type , 'agus' )
        pop_vals( : , 1 ) = 0.107;
        pop_vals_se( : , 1 ) = 0.001;
        pop_vals( xval >= 7.74 , 2 ) = 0.158;
        pop_vals_se( xval >= 7.74 , 2 ) = 0.002;
        pop_vals( xval < 6.77 , 3 )  = 0.099;
        pop_vals_se( xval < 6.77 , 3 )  = 0.003;
    elseif strcmp( type , 'bct' )
        pop_vals( : , 1 ) = 0;
        pop_vals_se( : , 1 ) = 0.001;
        pop_vals( xval >= 7.74 , 2 ) = 0.003;
        pop_vals_se( xval >= 7.74 , 2 ) = 0.002;
        pop_vals( xval < 6.77 , 3 ) = -0.008;
        pop_vals_se( xval < 6.77 , 3 ) = 0.004;
    elseif strcmp( type , 'bcg' )
        pop_vals( : , 1 ) = 0.018;
        pop_vals_se( : , 1 ) = 0.002;
        pop_vals( xval >= 7.74 , 2 ) = 0.048;
        pop_vals_se( xval >= 7.74 , 2 ) = 0.003;
        pop_vals( xval < 6.77 , 3 ) = -0.005;
        pop_vals_se( xval < 6.77 , 3 ) = 0.009;
    end
elseif Zid == 2
    pop_vals = NaN( neval , 4 );
    pop_vals_se = pop_vals;
    if strcmp( type , 'agus' )
        pop_vals( : , 1 ) = 0.107;
        pop_vals_se( : , 1 ) = 0.001;
        pop_vals( xval >= 3.5 , 2 ) = 0.097;
        pop_vals_se( xval >= 3.5 , 2 ) = 0.001;
        pop_vals( xval >= 2.5 & xval < 3.5 , 3 )  = 0.333;
        pop_vals_se( xval >= 2.5 & xval < 3.5 , 3 )  = 0.007;
        pop_vals( xval >= 2 & xval < 2.5 , 4 )  = 0.697;
        pop_vals_se( xval >= 2 & xval < 2.5 , 4 )  = 0.01;
    elseif strcmp( type , 'bct' )
        pop_vals( : , 1 ) = 0.000;
        pop_vals_se( : , 1 ) = 0.001;
        pop_vals( xval >= 3.5 , 2 ) = 0.000;
        pop_vals_se( xval >= 3.5 , 2 ) = 0.001;
        pop_vals( xval >= 2.5 & xval < 3.5 , 3 ) = -0.007;
        pop_vals_se( xval >= 2.5 & xval < 3.5 , 3 ) = 0.009;
        pop_vals( xval >= 2 & xval < 2.5 , 4 )  = 0.020;
        pop_vals_se( xval >= 2 & xval < 2.5 , 4 )  = 0.013;
    elseif strcmp( type , 'bcg' )
        pop_vals( : , 1 ) = 0.018;
        pop_vals_se( : , 1 ) = 0.002;
        pop_vals( xval >= 3.5 , 2 ) = 0.016;
        pop_vals_se( xval >= 3.5 , 2 ) = 0.002;
        pop_vals( xval >= 2.5 & xval < 3.5 , 3 )  = 0.028;
        pop_vals_se( xval >= 2.5 & xval < 3.5 , 3 )  = 0.016;
        pop_vals( xval >= 2 & xval < 2.5 , 4 )  = 0.179;
        pop_vals_se( xval >= 2 & xval < 2.5 , 4 )  = 0.023;
    end
end

%% Create figure
% Create graphs
f1 = figure('Position', get(0, 'Screensize'));
col = {'black','magenta','green','cyan'};
% First make the tau -> legends...
% RD effects
if Zid == 1
    plot( xv_t , tau , 'b','linewidth' , 2 );
    xlabel('School average transition score')
    leg_aux = {'Pop-Eleches et. al. - top tercile','Pop-Eleches et. al. - bottom tercile'};
    hold on
    % Pop-Elches
    for i = 1 : 3
       plot( xval , pop_vals( : , i ) , 'linewidth' , 2 , 'color', col{ i } , 'Linestyle', '-' ); 
    end
else
    errorbar( xv_t , tau , 2*tau_se, 'linewidth' , 2 , 'Linestyle', 'none','marker' , 'd' );
    xlabel('Number of schools in town')
    leg_aux = {'Pop-Eleches et. al. - 4 or more','Pop-Eleches et. al. - only 3','Pop-Eleches et. al. - only 2'};
    hold on
    % Pop-Elches
    for i = 1 : 4
        if i > 2
            errorbar( xval , pop_vals( : , i ) , 2 * pop_vals_se( : , i ) , 'linewidth' , 2 , 'color', col{ i } , 'Linestyle', '-' , 'marker' , 'o' ); 
        else
            errorbar( xval , pop_vals( : , i ) , 2 * pop_vals_se( : , i ) , 'linewidth' , 2 , 'color', col{ i } , 'Linestyle', '-' ); 
        end
    end
end

if Zid == 1
    % Create SE-s
    % For RD tree
    plot( xv_t , tau + 2*tau_se , '--', 'color','r','linewidth' , 2 );
    hold on
    plot( xv_t , tau - 2*tau_se , '--', 'color','r','linewidth' , 2 );
    for i = 1 : 3
        plot( xval , pop_vals( : , i ) + 2 * pop_vals_se( : , i ) , 'linewidth' , 2 , 'color', col{ i } , 'Linestyle', '--');
        plot( xval , pop_vals( : , i ) - 2 * pop_vals_se( : , i ) , 'linewidth' , 2 , 'color', col{ i } , 'Linestyle', '--' );
    end
end

ylabel('Treatment effect')
box('on');
grid('on');
legend('boxoff')
set(gca,'fontsize',22)
legends_cl = { 'RD tree','Pop-Eleches et. al. - ATE' , leg_aux{:} };
legend(legends_cl)

end
