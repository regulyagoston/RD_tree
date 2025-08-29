
function [ tau_p_mc , tse_p_mc , bias , hitCI , true_tau, bias_pctg ] = create_treatment_leaf( type , tL , tau_pred , tau_se_pred , ...
                                                          CIval )

MC = numel( tau_pred );

if tL == 1
    tau_p_mc = NaN( MC , 1 );
elseif tL == 2
    tau_p_mc = NaN( MC , 2 );
else
    tL = 5;
    tau_p_mc = NaN( MC , 5 );
end

tse_p_mc = tau_p_mc;
bias = tau_p_mc;
hitCI = tau_p_mc;

switch type
    case { 'binary-homogen-1tr' , 'binary-hetero-1tr' }
        true_tau = 1;
    case {  'binary-homogen-2tr' , 'binary-hetero-2tr' } 
        true_tau = [ 1 , -1 ];
    case { 'cnts-homogen-2tr' , 'cnts-hetero-2tr' } 
        true_tau = [ -2 , 2 ];
    case {'cnts-homogen-cntstr' , 'cnts-hetero-cntstr'}
        true_tau = [ -8 , -4 , 0 , 2 , 8 ];
    case 'rdrobust-M3-1tr'
        true_tau = 0.04;
    case 'rdrobust-M1-2tr'
        true_tau = [ 0.02 , 0.08 ];
    case 'rdrobust-M2-cntstr'
        true_tau = [ -5.45 , -6.45 , -7.45 , -8.45 , -9.45 ];
    otherwise
        error('No such DGP implemented')
end
true_tau = repmat( true_tau , [ MC , 1 ] );
        
for mc = 1 : MC
    if any( strcmp( type , { 'binary-homogen-1tr' , 'binary-hetero-1tr' ,'rdrobust-M3-1tr'} ) )
        tau_p_mc( mc ) = mean( tau_pred{ mc } );
        tse_p_mc( mc ) = mean( tau_se_pred{ mc } );
    elseif any( strcmp( type , { 'binary-homogen-2tr' , 'binary-hetero-2tr' } ) )
        tau_p_mc( mc , 1 ) = mean( tau_pred{ mc }( 1 : 2 ) );
        tau_p_mc( mc , 2 ) = mean( tau_pred{ mc }( 3 : 4 ) );
        tse_p_mc( mc , 1 ) = mean( tau_se_pred{ mc }( 1 : 2 ) );
        tse_p_mc( mc , 2 ) = mean( tau_se_pred{ mc }( 3 : 4 ) );
    elseif any( strcmp( type , { 'cnts-homogen-2tr' , 'cnts-hetero-2tr', 'rdrobust-M1-2tr'} ) )
        nPred = numel( tau_pred{ mc } ) ./ 2;
        nPred_se = max( [numel( tau_se_pred{ mc } ) ./ 2,1]);
        tau_p_mc( mc , 1 ) = mean( tau_pred{ mc }( 1 : nPred ) );
        tau_p_mc( mc , 2 ) = mean( tau_pred{ mc }( nPred+1 : end ) );        
        if nPred_se == 1
            tse_p_mc( mc , 1 ) = tau_se_pred{ mc };
            tse_p_mc( mc , 2 ) = tau_se_pred{ mc };
        else
            tse_p_mc( mc , 1 ) = mean( tau_se_pred{ mc }( 1 : nPred ) );
            tse_p_mc( mc , 2 ) = mean( tau_se_pred{ mc }( nPred+1 : end ) );
        end
        
    elseif any( strcmp( type , { 'cnts-homogen-cntstr' , 'cnts-hetero-cntstr','rdrobust-M2-cntstr' } ) )
        nPred = numel( tau_pred{ mc } ) ./ 5;
        nPred_se = max( [numel( tau_se_pred{ mc } ) ./ 2,1]);
        if nPred_se == 1
            for k = 1 : 5
                tau_p_mc( mc , k ) = mean( tau_pred{ mc }( ( k - 1 ) * nPred + 1 : k * nPred ) );
            end
            tse_p_mc( mc , : ) = tau_se_pred{ mc };
        else
            for k = 1 : 5
                tau_p_mc( mc , k ) = mean( tau_pred{ mc }( ( k - 1 ) * nPred + 1 : k * nPred ) );
                tse_p_mc( mc , k ) = mean( tau_se_pred{ mc }( ( k - 1 ) * nPred + 1 : k * nPred ) );
            end
        end
    else
        error('Not implemented')
    end
end

bias = true_tau - tau_p_mc;
bias_pctg = bias ./ tau_p_mc;

% Hit ratio for confidence intervals
lratio = ( 1 - CIval ) / 2;
CI_xs = norminv( [ lratio , CIval + lratio ] , 0 , 1 );

hitCI = tau_p_mc >= ( true_tau + CI_xs( 1 ) .* tse_p_mc ) & ...
        tau_p_mc <= ( true_tau + CI_xs( 2 ) .* tse_p_mc );

% sd_tau_p_mc = std( tau_p_mc );
% sd_tau_p_mc = mean( tse_p_mc );
% coef_se = repmat( coef_se' , [ 1 , tL ] ); 
% hitCI = tau_p_mc >= ( true_tau + CI_xs( 1 ) .* coef_se ) & ...
%         tau_p_mc <= ( true_tau + CI_xs( 2 ) .* coef_se );

end