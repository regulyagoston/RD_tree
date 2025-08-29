%% Cross validate to get optimal pruning/complexity parameter


function [ opt_beta , opt_oOS, s_oOS , beta , oOS ] = cross_validate( obj_tree , obj_sample , obj_optCART )


% 1st: Get the raw complexity parameters
[ prTrees , iOS , alpha ] = pruning( obj_tree , obj_sample , obj_optCART );

% 2nd: calculate 'typical values' instead of raw values
nM = numel( alpha );
beta = zeros( 1 , nM );
if any( alpha < 0 )
    %warning('There are negative candidate complexity parameters!')
    signA = sign( alpha );
else
    signA = ones( 1 , nM );
end
for i = 2 : nM - 1
    beta( i ) =  sqrt( abs( alpha( i ) .* alpha( i + 1 ) ) );
end
beta = beta .* signA;
beta( end ) = Inf;
if numel( unique( beta ) ) ~= numel( beta )
    error('cross_validate','Some problems with the complexity parameters! There are repetitions of the same value!')
end

%3rd: create cross-validation samples
[ train_sample , test_sample ] = create_CV_samples( obj_sample , obj_optCART , obj_optCART.numKfold );
oOS = NaN( obj_optCART.numKfold , nM );
e_tau_sq_1 = oOS;
e_tau_sq_2 = oOS;
bias_part = oOS;
var_part = oOS;
% Prune with betas and predict on test sample
prTrees_k   = repmat( tree , [ obj_optCART.numKfold , nM ] );
predTrees_k = repmat( tree , [ obj_optCART.numKfold , nM ] );
% % Original bandwidth for large tree
if isprop(obj_optCART,'bw' )
    orig_h = obj_optCART.bw;
    orig_b = obj_optCART.bw_b;
end
warnState = warning('off', 'all');  % Suppress all warnings temporarily
warning_flag = false(obj_optCART.numKfold,nM);
    

% try
%     if ~obj_optCART.paralell
%         error('This is only a shortcut if paralell processing is not allowed!')
%     end
%     parfor k = 1 : obj_optCART.numKfold
%         % If there is bandwidths modify them to fit train sample
%         if isprop(obj_optCART.bw )
%             h = modify_bandwidth( obj_optCART.bwselect, obj_optCART.bw, ...
%                 obj_sample.n_tr, train_sample( k ).n_tr , obj_optCART.orderPolinomial );
%             b = modify_bandwidth( obj_optCART.bwselect, obj_optCART.bw_b, ...
%                 obj_sample.n_tr, train_sample( k ).n_tr , obj_optCART.orderPolinomial );
%             obj_optCART = copyWithChange(obj_optCART,'bw', h);
%             obj_optCART = copyWithChange(obj_optCART,'bw_b', b);
%         end
%         % Fit a full model on the cv data set
%         obj_tree_k = growRDDtree( train_sample( k ) , obj_optCART );
%         % For each candidate complexity parameter
%         for i = 1 : nM
%             % Prune trees based on training sample
%             prTrees_k( k , i )   = pruning( obj_tree_k , train_sample( k ) , obj_optCART , beta( i ) );
%             % Use the pruned trees to predict or estimate on test sample
%             predTrees_k( k , i ) = clearTree( prTrees_k( k , i ) );
%             % If there is bandwidths modify them to fit estimation sample
%             if isprop(obj_optCART.bw )
%                 h = modify_bandwidth( obj_optCART.bwselect, obj_optCART.bw, ...
%                     obj_sample.n_tr, test_sample( k ).n_est , obj_optCART.orderPolinomial );
%                 b = modify_bandwidth( obj_optCART.bwselect, obj_optCART.bw_b, ...
%                     obj_sample.n_tr, test_sample( k ).n_est , obj_optCART.orderPolinomial );
%                 obj_optCART = copyWithChange(obj_optCART,'bw', h);
%                 obj_optCART = copyWithChange(obj_optCART,'bw_b', b);
%             end
%             if any( strcmp( obj_optCART.type , 'prediction' ) )
%                 predTrees_k( k , i ) = predictTree( predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k ) , obj_optCART );
%             elseif any( strcmp( obj_optCART.type , {'CATE','CLATE'} ) )
%                 predTrees_k( k , i ) = updateTree(  predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k )   , obj_optCART , true );
%             else
%                 error('cross_validate:wrongInputType',['There are no such type implemeted as: ' obj_optCART.type ] )
%             end
%             % Calcilate out of sample criterion on test sample
%             p_l_ID        = findleaves( predTrees_k( k , i ) );
%             oOS( k , i )  = obj_optCART.crit_oos( [ predTrees_k( k , i ).Node{ p_l_ID } ] , test_sample( k ) , obj_optCART );
%         end
%     end
% catch
    for k = 1 : obj_optCART.numKfold
        % Clear any previous warning
        lastwarn('');
        % % If there is bandwidths modify them to fit train sample
        if isprop(obj_optCART,'bw' )
            h = modify_bandwidth( obj_optCART.bwselect, orig_h, ...
                obj_sample.n_tr, train_sample( k ).n_tr , obj_optCART.orderPolinomial );
            b = modify_bandwidth( obj_optCART.bwselect, orig_b, ...
                obj_sample.n_tr, train_sample( k ).n_tr , obj_optCART.orderPolinomial );
            obj_optCART = copyWithChange(obj_optCART,'bw', h);
            obj_optCART = copyWithChange(obj_optCART,'bw_b', b);
        end
        % Fit a full model on the cv data set
        obj_tree_k = growRDDtree( train_sample( k ) , obj_optCART );
        % For each candidate complexity parameter
        for i = 1 : nM
            % Prune trees based on training sample
            prTrees_k( k , i )   = pruning( obj_tree_k , train_sample( k ) , obj_optCART , beta( i ) );
            % Use the pruned trees to predict or estimate on test sample
            predTrees_k( k , i ) = clearTree( prTrees_k( k , i ) );
            % % If there is bandwidths modify them to fit estimation sample
            if isprop(obj_optCART,'bw' )
                h = modify_bandwidth( obj_optCART.bwselect, orig_h, ...
                    obj_sample.n_tr, test_sample( k ).n_est , obj_optCART.orderPolinomial );
                b = modify_bandwidth( obj_optCART.bwselect, orig_b, ...
                    obj_sample.n_tr, test_sample( k ).n_est , obj_optCART.orderPolinomial );
                obj_optCART = copyWithChange(obj_optCART,'bw', h);
                obj_optCART = copyWithChange(obj_optCART,'bw_b', b);
            end
            if any( strcmp( obj_optCART.type , 'prediction' ) )
                predTrees_k( k , i ) = predictTree( predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k ) , obj_optCART );
            elseif any( strcmp( obj_optCART.type , {'CATE','CLATE'} ) )
                predTrees_k( k , i ) = updateTree(  predTrees_k( k , i ) , prTrees_k( k , i ) , test_sample( k ) , obj_optCART , true );
            else
                error('cross_validate:wrongInputType',['There are no such type implemeted as: ' obj_optCART.type ] )
            end
            % Calcilate out of sample criterion on test sample
            p_l_ID        = findleaves( predTrees_k( k , i ) );
            [ oOS( k , i ), e_tau_sq_1( k, i ), e_tau_sq_2( k, i ), ...
                bias_part( k, i ), var_part( k, i ) ]  = obj_optCART.crit_oos( [ predTrees_k( k , i ).Node{ p_l_ID } ] , test_sample( k ) , obj_optCART );
            %
        end
        % --- Check for warning ---
        [msg, ~] = lastwarn;
        if ~isempty(msg)
            warning_flag( k, i ) = true;
        end
    end
%end

% If there were a warning neglect that fold
oOS( warning_flag ) = NaN;

warning(warnState);                 % Restore original warning state
% if ( any( warning_flag ) )
%     warning('During cross-validation, there were instances where warning messages had been suppressed and the accompanying out-of-sample criterion values neglected!')
% end

if obj_optCART.cv_rm_outlier
    %for k = 1:nM
        log_outlier = isoutlier( oOS(:), 'gesd' );%(:,k) ); , 'percentiles', [0.001, 0.999 ]
        %oOS( log_outlier,k ) = NaN;
        log_outlier = reshape( log_outlier, [ obj_optCART.numKfold, nM ] );
        oOS( log_outlier ) = NaN;
    %end    
    s_oOS = mean(oOS,1, 'omitnan');   
    if any( sum( log_outlier, 1 ) == obj_optCART.numKfold )
        s_oOS( sum( log_outlier, 1 ) == obj_optCART.numKfold ) = Inf;
    end
else
    % Get the mean of out-of-sample values
    s_oOS = sum( oOS, 'omitnan' );
end


% Get the minimum value
[ mVal , mID ] = min( s_oOS );
% If 1SE rule
if obj_optCART.CV1SE
    u_range  = mVal + ( std( oOS( : , mID ), 'omitnan' ) ./ sum( ~isnan( oOS( : , mID ) )) );
    val_id   = s_oOS <= u_range;
    aux_beta = beta( val_id );
    opt_beta = aux_beta( end );
else
    opt_beta = beta( mID );
end
opt_oOS = mVal;


end