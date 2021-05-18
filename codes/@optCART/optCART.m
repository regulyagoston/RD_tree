classdef optCART < dynamicprops
    %% Growing a large tree
    properties
        % Maximum number of nodes to create a large tree
        maxNodes = Inf;
        % Maximal level of deepth of the tree
        maxLevel = Inf;
        % Maximal number of leaves
        maxLeaves = 50;
        % Minimal observation that a leaf needs to contain
        minObs = 200;
        % Complexity parameter for growing a tree -> if IOS does not
        % increase by this amount there is no further splitting. Only valid
        % if overfit_by_crit is false
        cp = 0.001;
        % Maximum iteration while growing a tree
        maxIterGrow = 200;
        % Set seed for the random number generation
        seed = 777;
    end
    %% Splitting options to find the cut value, when splitting the covariates
    properties 
        % Splitting type
        splitType = 'causal-trees-RDD'; 
                       % 'pctg_wM' - using percentiles with median forced
                       % 'pctg' - using percentiles without forcing median
                       % 'uni_support' - using uniformly spaced values
                       % 'rpart' - using all the possible splitting value (same as Hastie et. al.)
                       % 'causal-trees-AI'  - Athey-Imbens, using buckets
                       % 'causal-trees-RDD' - splitting optimized for RDD, using buckets
        % Maximum split that the algorithm computes by hand
        numSplit  = 10;
    end 
    %% Criterion type and order of polinomials used for estimation
    properties
       type       = 'CATE';      % Type can be:
                                      %     'prediction'   - focuses on predictive problems
                                      %     'CATE'         - estimates the conditional average treatment effect
                                      %     'CLATE'        - estimates the conditional local average treatment effect
                                      
       model      = 'mean'            % Model can be:
                                      %     'mean'   - estimate simple means of the outcome
                                      %     'linear' - estimate linear regression with OLS based on given X covariates
                                      
       eval_level = 'leaf-by-leaf'    % Evaluation level can be:
                                      %     'leaf-by-leaf' - everything is calculated on a node level:
                                      %                         once it is calculated it is saved
                                      %     'tree'         - evaluation is based on the whole tree:
                                      %                         it calculates everything each time on the whole tree    
                                      
       criteria   = 'RDD';            % Criterion can be:
                                      %     'athey-imbens' - uses EMSE criterion for experimental/unconfounded treatment as in Athey-Imbens (2006)
                                      %     'RDD'          - uses EMSE criterion for RDD treatment developed in Reguly (2020)
                                      %     'MSE'          - uses the MSE for predicted values
                                      
       honest     = true;             % Honest tree or simple additive tree
                                      %     use honest evaluation if ture 
                                      %     if false it uses the classical adaptive tree   
    end
    %% Self-functions
    properties %( Access = private )
       crit_split   % Gives the splitting criterion function - used to find the best split value for each feature
       crit_ios     % Gives the in-sample criterion function - used to prune the tree
       crit_oos     % Gives the splitting criterion function - used to find complexity parameter
       split_fnc    % Gives the splitting function           - gives the split points
    end
    %% Other properties
    properties
        % Number of folds during cross validation
        numKfold = 5;
        % Scale the features
        scaleFeatures = false;
        % Cross-validation use 1-SE rule to get optimal pruning parameter
        CV1SE = true;
        % Use paralell processing
        paralell = false;
    end
    
    methods
        function obj = optCART()
            return;
        end
        
        function set.honest( obj , val )
           if ~islogical( val )
               error('Invalid value for honest, use logical: true is to use honest approach')
           else
               obj.honest = val;
           end
           % Add dynamic properties in case of honest estimation
           if obj.honest
               % Add Leave-one-out estimator possibility if honest estimation
               p0 = addprop( obj , 'LOO' );
               p0.SetMethod = @set_LOO;
               obj.LOO = false;
               % Criteria weight -> if set to 1 it neglects the variance part
               p1 = addprop( obj , 'criteria_weight' );
               p1.SetMethod = @set_criteria_weight;
               obj.criteria_weight = 0.5;
           end
           
        end
        function set.type( obj , val )
            if ~any( strcmp( val , { 'prediction','CATE','CLATE' } ) ) 
                error('optCART:type:incorrectType:val',...
                    ['Invalid value for type of model, use one of the following: \n',...
                     'CATE   - Conditional Average Treatment effect estimation \n',...
                     'CLATE  - Conditional Local Average Treatment effect estimation - uses first stage EQ for getting treatment \n',...
                     'prediction - predict the outcome variable, no treatment effect'])
            else
                obj.type = val;
                if strcmp( obj.type , 'CLATE' )
                   p1 = addprop( obj , 'CLATE_crit_tau_w' );
                   p1.SetMethod = @set_CLATE_crit_tau_w;
                   obj.CLATE_crit_tau_w = 0.0001;
               end
            end
        end
        function set.model( obj , val )
            if ~any( strcmp( val , { 'mean','linear' } ) ) 
                error('optCART:model:incorrectType:val',...
                    ['Invalid value for model, use one of the following: \n',...
                    'mean   - heterogeneous treatment parameter estimation \n',...
                    'linear - estimate linear regression with OLS based on given X covariates'])
            else
                obj.model = val;
            end
            % Add dynamic properties in case of linear regression model
            if strcmp( obj.model , 'linear' )
               p0                 = addprop( obj , 'orderPolinomial' );
               p0.SetMethod       = @set_orderPolinomial;
               obj.orderPolinomial = 1;
               p1                 = addprop( obj , 'varEstimator' );
               p1.SetMethod       = @set_varEstimator;
               obj.varEstimator   = 'hce-1';
            end
        end
        function set.eval_level( obj , val )
            if ~any( strcmp( val , { 'leaf-by-leaf','tree' } ) ) 
                error('optCART:eval_level:incorrectType:val',...
                    ['Invalid value for eval_level, use one of the following: \n',...
                    'leaf-by-leaf  - everything is calculated on a node level \n',...
                    'tree          - evaluation is based on the whole tree'])
            else
                obj.eval_level = val;
            end
        end
        function set.criteria( obj , val )
            if ~any( strcmp( val , { 'MSE','athey-imbens','RDD' } ) ) 
                error('optCART:criteria:incorrectType:val',...
                    ['Invalid value for criteria, use one of the following: \n',...
                    'athey-imbens - treatment EMSE criterion based on athey-imbens \n',...
                    'RDD          - treatment EMSE criterion for RDD \n',...
                    'MSE          - errors from model fit is used to calculate MSE values'])
            else
                obj.criteria = val;
            end
            if strcmp( val , 'RDD' ) 
               p0 = addprop( obj , 'scale_runningvariable_w_c_2_zero' );
               p0.SetMethod = @set_scale_runningvariable_w_c_2_zero;
               obj.scale_runningvariable_w_c_2_zero = true;
            end
        end
        function set.splitType( obj , val )
            if ~any( strcmpi( val , {'pctg_wM','pctg','uni_support','rpart','causal-tree-AI','causal-tree-RDD'} ) )
                error('optCart:splitType:wrongInput',...
                    ['Invalid value for split type, use one of the following value: \n',...
                    'pctg_wM         - using percentiles with median forced \n',...
                    'pctg            - using percentiles without forcing median \n',...
                    'uni_support     - using uniformly spaced values between minimum and maximum values \n',...
                    'rpart           - using all the possible splitting value (same as Hastie et. al.) \n',...
                    'causal-tree-AI  - algorithm develpoled based on Athey-Imbens, using buckets \n' , ...
                    'causal-tree-RDD - algorithm develpoled based on Reguly, using buckets']);
            else
                % Set the splitting type and the function
                obj.splitType = val;
                setSplitFnc( obj );
            end
            % Add dynamic properties in case of casual-trees splitting criterion
            if strncmpi( val , 'causal-tree' , 11 )
                % number of observations in a bucket enforced by splitting
               p0 = addprop( obj , 'obsBucket' );
               p0.SetMethod = @set_obsBucket;
               obj.obsBucket = 5;
               % maximum number of bucket to use
               p1 = addprop( obj , 'maxBucket' );
               p1.SetMethod = @set_maxBucket;
               obj.maxBucket = 100;
            end
        end
        function set.numKfold( obj , val )
           if numel( val ) ~= 1 || val < 0 || mod( val , 1 ) ~= 0
               error('optCart:numKfold:wrongInput',...
                     'Invalid value for number of K-fold, use a positive integer, greater than 1')
           else
               obj.numKfold = val;
           end
        end        
        function set.maxNodes( obj , val )
            if numel( val ) ~= 1 || val < 0 || ( ~isinf( val ) && mod( val , 1 ) ~= 0 )
               error('optCart:maxNodes:wrongInput',...
                     'Invalid value for maximum number of nodes, use an integer larger than 0 or inf.')
           else
               obj.maxNodes = val;
           end
        end
        function set.maxLevel( obj , val )
            if numel( val ) ~= 1 || val < 0 || ( ~isinf( val ) && mod( val , 1 ) ~= 0 )
               error('optCart:maxLevel:wrongInput',...
                     'Invalid value for maximum level reached by the tree, use an integer equal or larger than 0')
           else
               obj.maxLevel = val;
           end
        end
        function set.maxLeaves( obj , val )
            if numel( val ) ~= 1 || val < 0 || ( ~isinf( val ) && mod( val , 1 ) ~= 0 )
               error('optCart:maxLeaves:wrongInput',...
                     'Invalid value for maximum leaves grown by the tree, use an integer larger than 0')
           else
               obj.maxLeaves = val;
           end
        end
        function set.minObs( obj , val )
            if numel( val ) ~= 1 || val < 0 || mod( val , 1 ) ~= 0
               error('optCart:minObs:wrongInput',...
                     'Invalid value for minimum number of observations in each leaf, use an integer larger than 0')
           else
               obj.minObs = val;
           end
        end
        function set.maxIterGrow( obj , val )
            if numel( val ) ~= 1 || val < 0 || ( ~isinf( val ) && mod( val , 1 ) ~= 0 )
               error('optCart:minObs:wrongInput',...
                     'Invalid value for maximum number of iteration used during growing the tree, use an integer larger than 0')
            else
               if isinf( val )
                   val = 10^6;
               end
               obj.maxIterGrow = val;
           end
        end
        function set.numSplit( obj , val )
            if numel( val ) ~= 1 || val < 0 || mod( val , 1 ) ~= 0
               error('optCart:numSplit:wrongInput',...
                     'Invalid value for maximum number of splits during searching for the splits, use an integer larger than 0')
           else
               obj.numSplit = val;
           end
        end
        function set.cp( obj , val )
           if numel( val ) ~= 1 || val < 0 || isinf( val )
               error('optCart:cp:wrongInput',...
                     [ 'Invalid value for complexity parameter, use a positive pr 0 value!' ] )
           else
               obj.cp = val;
           end
        end  
        function set.scaleFeatures( obj , val )
           if ~islogical( val )
               error('optCart:scaleFeatures:wrongInput',...
                     ['Invalid value for scaleFeatures, use logical: \n',...
                      'true  - standardize each feature to have zero mean and standard deviation 1 \n',...
                      'false - no scaling for features.'])
           else
               obj.scaleFeatures = val;
           end
        end
        function set.CV1SE( obj , val )
           if ~islogical( val )
               error('optCart:CV1SE:wrongInput',...
                     ['Invalid value for cross-validation 1-SE rule (CV1SE), use logical: \n',...
                      'true  - use 1-SE rule to get optimal pruning parameter  \n',...
                      'false - parameter is given by the minimum value.'])
           else
               obj.CV1SE = val;
           end
        end
        function set.seed( obj , val )
            if numel( val ) ~= 1 || val < 0 || mod( val , 1 ) ~= 0
               error('optCart:seed:wrongInput',...
                     'Invalid value for seed of random number generation, use an integer larger than 0')
           else
               obj.seed = val;
            end
        end
        function set.paralell( obj , val )
           if ~islogical( val )
               error('Invalid value for paralell, use logical: true is to use paralell processing during the esimation!')
           else
               obj.paralell = val;
           end
        end
    end
    
    %% Dynamic properties set functions
    methods
        % Leave-one-out estimator
        function set_LOO( obj , val )
            if ~islogical( val ) || ( val && ~obj.honest ) || numel( val ) ~= 1
                error('optCart:set_LOO:wrongInput',...
                      'Invalid value for LOO, use logical: true is to use leave-one-out estimator. Only working with honest approach!')
            else
                obj.LOO = val;
            end
        end
        % Criteria weight used in evaluation
        function set_criteria_weight( obj , val )
           if ~isnumeric( val ) || numel( val ) ~= 1 || val < 0 || val > 1
               error('optCart:set_criteria_weight:wrongInput',...
                     'Invalid value for criteria_weight, use scalar between 0-1')
           else
               obj.criteria_weight = val;
           end
        end
        % Number of observations in a bucket
        function set_obsBucket( obj , val )
           if numel( val ) ~= 1 || val < 0 || mod( val , 1 ) ~= 0
               error('optCart:set_obsBucket:wrongInput',...
                     'Invalid value for number of observations in a bucket, use an integer larger than 0')
           else
               obj.obsBucket = val;
           end
        end
        % Maximum number of buckets used
        function set_maxBucket( obj , val )
           if numel( val ) ~= 1 || val < 0 || ( ~isinf( val ) && mod( val , 1 ) ~= 0 )
               error('optCart:set_maxBucket:wrongInput',...
                     'Invalid value for maximum number of buckets, use an integer larger than 0')
           else
               obj.maxBucket = val;
           end
        end
        % Set the order of polinomial
        function set_orderPolinomial( obj , val )
            if ~( numel( val ) == 1 && mod( val , 1 ) == 0 && val > 0 && val < 7 )
                error('optCart:orderPolinomial:wrongInput',...
                      'Invalid value order of polinomial for X, use a positive integer, less than 7!');
            else
                % Set the order of polinomial
                obj.orderPolinomial = val;
            end            
        end
        % Set the variance Estimator
        function set_varEstimator( obj , val )
            if ~any( strcmpi( val , {'simple','hce-0','hce-1'} ) )
                error('optCart:varEstimator:wrongInput',...
                     ['Invalid value for variance estimator! Choose one of the following: \n',...
                      'simple: simple homoskedastic estimator \n',...
                      'hce-0: White (1980) heteroscedasticity-consistent estimator \n',...
                      'hce-1: hce-0 with degree of freedom correction (Hinkley, 1977) \n']);
            else
                % Set the order of polinomial
                obj.varEstimator = val;
            end            
        end
        % Use scaling for the running variable
        function set_scale_runningvariable_w_c_2_zero( obj , val )
            if ~islogical( val ) || numel( val ) ~= 1
                error('optCart:scale_runningvariable_w_c_2_zero:wrongInput',...
                      [ 'Invalid value for scale_runningvariable_w_c_2_zero, use a logical: \n',...
                        'true  - scale running variable such that c=0 will be the new cut-off value\n',...
                        'false - no scaling'])
            else
                obj.scale_runningvariable_w_c_2_zero = val;
            end
        end
        % Set the critical value for the first stage in CLATE
        function set_CLATE_crit_tau_w( obj , val )
            if ~( numel( val ) == 1 && val > 0 && val < 1 )
                error('optCart:CLATE_crit_tau_w:wrongInput',...
                      'Invalid value for critical value for CLATE restriction (first stage effect), use an integer between 0 and 1!');
            else
                % Set the critical value for the first stage in CLATE
                obj.CLATE_crit_tau_w = val;
            end            
        end
    end
    
    methods
        setSplitFnc( obj );
        setCritFnc( obj );
    end
end

