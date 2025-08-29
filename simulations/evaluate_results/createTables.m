%% Create pretty tables
clear all
clc

types = {'rdrobust-M3-1tr','rdrobust-M1-2tr','rdrobust-M2-cntstr'};
algo = {'np_uni_','forest_uni_'};
path = ''

MC = 1000;
ci_val = 0.95;
fuzzy = true;

isSEused = false; % Do not change (if 1SE rule usef or not)
lBW = false; % Do not change (parametric vs nonparametric comparison)

if strcmp( types{1}(1:8) , 'rdrobust' )
    sigmaY0 = 0.05^2;%[0.05,0.1295].^2;
    N  = [1000;5000;10000;50000];
else
    sigmaY0 = 1;%[ 0.1 , 1, 10 ];
    N  = [1000,5000,10000,50000];
end

sigmaY = sigmaY0;
nT = numel( types );
if numel( sigmaY ) > 1
    nS = numel( sigmaY );
    header_ch = '$\sigma_\epsilon$';
else
    nS = numel( N );
    header_ch = 'Num obs.';
    if sigmaY < 1
        aux_s2 = num2str( sigmaY );
        sStr = horzcat('0',aux_s2( 3 : end ) );
    else
        sStr = num2str( sigmaY );
    end
end

if fuzzy
    id_tab = 'fRDD_';
else
    id_tab = 'sRDD_';
end
if lBW
    id_bw = '_BW';
else
    id_bw = '';
end



%% Creates Table with InfMSE, number of leaves, DGP found
charTable = '';%horzcat( '\begin{tabular}{|l|l|c|c|} \hline DGP & ' , header_ch , ' & inf. MSE & $\# \Pi$ \hline' );
for  i = 1 : nT
    charTable = horzcat( charTable , '\multirow{3}{*}{DGP-' , num2str( i ) , '} ' );
    for j = 1 : nS
        if numel( sigmaY ) > 1
            if sigmaY( j ) < 1
                aux_s2 = num2str( sigmaY( j ) );
                sStr = horzcat('0',aux_s2( 3 : end ) );
            else
                sStr = num2str( sigmaY( j ) );
            end
            N_ch = num2str( N );
        else
            N_ch = num2str( N( j ) );
        end
        
        % Select the proper file
        if isSEused
            cv_str = '_1SE';
        else
            cv_str = '_noSE';%_rho1_no_other';%_arch';
        end
        for a = 1:numel(algo)
            try
                % Load the file
                filename = horzcat( path , cv_str(2:end) , '/' , id_tab ,algo{ a }, types{ i } , '_' , sStr ,'_' , N_ch , cv_str,
                     id_bw );
                load( filename );
                sigmaY = sigmaY0;
                
                
                
                % Estimate statistics for treatment effects
                tL = str2double( types{ i }( end - 2 ) );
                [ ~ , ~ , bias_j , hitCI , true_tau, bias_pctg ] = create_treatment_leaf( ...
                    types{ i } , tL , tau_pred(1:MC) , tau_se_pred(1:MC) , ci_val );
                nTau = size( bias_j , 2 );
                
                % Set the names for tables
                if numel( sigmaY ) > 1
                    if sigmaY( j ) == 1
                        sigStr = '$   & ';
                    elseif sigmaY( j ) == 10
                        sigStr = '$  & ';
                    else
                        sigStr = '$ & ';
                    end
                    colVar = horzcat( '& $ \sigma = ' , num2str( sigmaY( j ) ) );
                else
                    sigStr = '$ & ';
                    colVar = horzcat( '& $ N = ' , num2str( N( j ) ) );
                end
                
                if numel( algo{a} ) >= 7 & strcmp( algo{a}(1:7), 'forest_')
                    % Creating table for performance
                    charTable = horzcat( charTable ,  ...
                        num2str( mean( inf_mse(1:MC) ) , '%1.4f' ) , ' & ', ...
                        num2str( mean( bias(1:MC) ) , '%1.4f' ) , ' & ', ...
                        num2str( mean( coverage(1:MC) ) , '%1.4f' ) , ' \\ ' );
                else
                    % Creating table for performance
                    charTable = horzcat( charTable , colVar , sigStr  , ...
                        num2str( mean( inf_mse(1:MC) ) , '%1.4f' ) , ' & ', ...
                        num2str( mean( bias(1:MC) ) , '%1.4f' ) , ' & ', ...
                        num2str( mean( coverage(1:MC) ) , '%1.4f' ) , ' & ', ...
                        num2str( mean( numLeaf(1:MC) ) , '%1.2f' ) , ' & ' );
                end      
            catch
                % Creating table for performance
                    charTable = horzcat( charTable ,  ...
                        num2str( NaN , '%1.4f' ) , ' & ', ...
                        num2str( NaN , '%1.4f' ) , ' & ', ...
                        num2str( NaN , '%1.4f' ) , ' \\ ' );
            end
        end
        % Add last line
        if  j == nS
            charTable = horzcat( charTable , ' \hline' );
        end
    end
end

charTable
