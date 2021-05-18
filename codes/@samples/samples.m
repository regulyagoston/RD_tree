classdef samples < dynamicprops
   
   %% Sample observations
   properties ( SetAccess = protected )
      % Training samples
      Y_tr  % Outcome variable
      Z_tr  % Features for tree building
   end
   
   %% Dependent variables of the class
   properties ( Dependent )
       n_tr      % Number of observation in training sample
       n_est     % Number of observation in estimation sample (conditional)
       K         % Number of features in Z
       L         % Number of regressors in X (conditional)
   end
   properties
      train_to_est_ratio = 1; 
   end
   properties ( Hidden )
      log_X_tr   = false;
      log_W_tr   = false;
      log_ITT_tr = false;
      log_cl_tr  = false;
   end
   
   
    methods
        function obj = samples( obj_optCART , Y , Z , varargin )
            
            %% Check for basic inputs
            if ~isa( obj_optCART , 'optCART' )
                error('samples:WrongInput:obOptCART',...
                      'objOptCART must be an optCART object!')
            end
            if all( Y( 1 ) == Y )
                error('samples:WrongInput:Y' , ...
                       'Invalid input there is no variation in Y use a random variable!' )
            end
            if any( all( Z( 1 , : ) == Z ) )
                Z( : , all( Z( 1 , : ) == Z ) ) = [];
                warning('samples:WrongInput:Z' , ...
                        'There is no variation in one or more variables of the Zs, these variables are dropped!' )
            end
            n = numel( Y );
            if n ~= size( Z , 1 )
                error('samples:WrongInput' , ...
                      'Number of observations in Y and Z differs from each other!' )
            end
            if any( isnan( Y ) )
                error('samples:WrongInput:Y' , ...
                      'There are NaN values in Y! First clean your data!' )
            end
            if any( any( isnan( Z ) ) )
                error('samples:WrongInput:Z' , ...
                      'There are NaN values in Z! First clean your data!' )
            end
            %% Input-parser
            iP = inputParser;
            iP.CaseSensitive = true;
            if obj_optCART.honest
                addParameter( iP , 'index4TrEst' , [] );
            end
            if any( strcmp( obj_optCART.type , {'CATE','CLATE'} ) )
                addParameter( iP , 'W' , [] );
            end
            if strcmp( obj_optCART.model , 'linear' )
                addParameter( iP , 'X' , [] );
                addParameter( iP , 'cluster' , [] );
            end
            if strcmp( obj_optCART.criteria , 'RDD' )
                addParameter( iP , 'cutoff' , [] );
            end
            parse( iP , varargin{:} );
            
            %% Check for the inputParser's inputs
            % Indexes for training and estimation sample
            if obj_optCART.honest
                ind_TrEst = iP.Results.index4TrEst;
                if isempty( ind_TrEst ) 
                    error('samples:WrongInput:index4TrEst' , ...
                          ['No values were given for indexes for training and estimation sample! \n',...
                          'In case of honest estimation you must define such index matrix!'])
                elseif ~all( unique( ind_TrEst( : ) )' == 1 : n )
                    error('samples:WrongInput:index4TrEst' , ...
                        [ 'Invalid input for indexes for training and estimation! Use all integer values for the indexes between 1 and '  num2str( n ) '!' ] )
                elseif size( ind_TrEst , 2 ) ~= 2
                    error('samples:WrongInput:index4TrEst' , ...
                          'Invalid input for indexes for training and estimation! Use: n times 2 matrix.' )
                end
                ind_tr  = unique( ind_TrEst( : , 1 ) , 'stable' );
                ind_est = unique( ind_TrEst( : , 2 ) , 'stable' );                
                obj.Y_tr  = Y( ind_tr , : );
                obj.Z_tr  = Z( ind_tr , : );
                % Add protected properties
                p1 = addprop( obj , 'Y_est' );
                p2 = addprop( obj , 'Z_est' );
                p1.SetAccess = 'protected';
                p2.SetAccess = 'protected';
                obj.Y_est = Y( ind_est , : );
                obj.Z_est = Z( ind_est , : );
                % train to estimation ratio
                %p3 = addprop( obj , 'train_to_est_ratio' );
                %p1.SetAccess = 'protected';
                obj.train_to_est_ratio = numel( obj.Y_tr ) ./ numel( obj.Y_est );
            else
                obj.Y_tr  = Y;
                obj.Z_tr  = Z;
            end
            
            % Add treatment variable if not sharp RDD
            if any( strcmp( obj_optCART.type , { 'CATE' , 'CLATE' } ) ) && ...
                ~( strcmp( obj_optCART.type , 'CATE' ) && strcmp( obj_optCART.criteria , 'RDD' ) )
            
                W = iP.Results.W;
                if isempty( W ) 
                     error( 'samples:WrongInput:W' ,...
                           'No treatment variable is given! For treatment type CART you need to define a treatment variable W!' )
                elseif ~islogical( W ) || ~all( size( W ) == size( Y ) )
                    error( 'samples:WrongInput:W' ,...
                           'Invalid type and size for treatment, use same sample size as Y and logicals!' )
                end
                if any( isnan( W ) )
                    error('samples:WrongInput:W' , ...
                          'There are NaN values in W! First clean your data!' )
                end
                % Variables for estimation
                p1 = addprop( obj , 'W_tr' );
                p1.SetAccess = 'protected';
                p2 = addprop( obj , 'p_tr' );
                p2.SetAccess = 'protected';
                p3 = addprop( obj , 'ITT_tr' );
                p3.SetAccess = 'protected';
                if obj_optCART.honest
                    p4 = addprop( obj , 'W_est' );
                    p4.SetAccess = 'protected';
                    p5 = addprop( obj , 'p_est' );
                    p5.SetAccess = 'protected';
                    p6 = addprop( obj , 'ITT_est' );
                    p6.SetAccess = 'protected';
                    obj.W_tr  = W( ind_tr , : );
                    obj.W_est = W( ind_est , : );
                    obj.p_est = sum( obj.W_est ) ./ obj.n_est;
                    obj.ITT_tr  = W( ind_tr , : );
                    obj.ITT_est = W( ind_est , : );
                else
                    obj.W_tr    = W;
                    obj.ITT_tr  = W;
                end
                % Note: ITT: intent-to-treat and in these cases these are
                % the same
                obj.p_tr  = sum( obj.W_tr ) ./ obj.n_tr;
            end
            
            if strcmp( obj_optCART.model , 'linear' ) || strcmp( obj_optCART.criteria , 'RDD' )
                X = iP.Results.X;
                if isempty( X ) 
                    error('samples:WrongInput:X' , ...
                           'No variable for X is given! For linear model or RDD criterion you need to define X variable(s)!' )
                end
                if any( all( X( 1 , : ) == X ) )
                    error('samples:WrongInput:X' , ...
                           'Invalid input there is no variation in X, use a random variable!' )
                end
                if size( X , 1 ) ~= n
                    error('samples:WrongInput:X' , ...
                          [ 'Invalid input for X, L times ' , num2str( n ) , ' matrix!' ] )
                end
                if strcmp( obj_optCART.criteria , 'RDD' )
                    if size( X , 2 ) ~= 1
                        warning('samples:WrongInput:X' , ...
                              'In case of RDD you must provide your running variable as the first vector in X!' )
                    end
                end
                if any( isnan( X ) )
                    error('samples:WrongInput:X' , ...
                          'There are NaN values in X! First clean your data!' )
                end
                if obj_optCART.orderPolinomial ~= 1
                    % Only takes polinomials for the first vector
                    X1p = createPolinomials( X( : , 1 ) , obj_optCART.orderPolinomial );
                    X = [ X1p , X( : , 2 : end ) ];
                end
                p1 = addprop( obj , 'X_tr' );
                p1.SetAccess = 'protected';
                if obj_optCART.honest
                    p2 = addprop( obj , 'X_est' );
                    p2.SetAccess = 'protected';
                    obj.X_tr  = X( ind_tr , : );
                    obj.X_est = X( ind_est , : );
                else
                    obj.X_tr  = X;
                end
                % Set the number of covariates parameter
                p3 = addprop( obj , 'nCov' );
                p3.SetAccess = 'protected';
                obj.nCov = size( obj.X_tr , 2 );
            end
            
            % In case of RDD
            if strcmp( obj_optCART.criteria , 'RDD' )
                % Cut-off value
                c = iP.Results.cutoff;
                if isempty( c )
                    error('samples:WrongInput:c' , ...
                          'No value for cut-off parameter is given! Define the cutoff scalar for RDD!' )
                elseif ~isnumeric( c ) || numel( c ) ~= 1 || isinf( c ) || c < min( X( : , 1 ) ) || c > max( X( : , 1 ) )
                    error('samples:WrongInput:c' , ...
                          [ 'Invalid input for cutoff value, use a scalar \n',...
                           'smaller than max value of the first vector in X: ' num2str( max( X( : , 1 ) ) ) ,'\n' , ...
                           'and larger value than minimum of the first vector in X: ' , num2str( min( X( : , 1 ) ) ) ] );
                end
                p = addprop( obj , 'c' );
                p.SetAccess = 'protected';
                obj.c  = c;
                % Set intent-to-treatment variable -> always use the first given vector!    
                ITT = X( : , 1 ) >= c;
                % In case of sharp RDD create and set the variable
                if strcmp( obj_optCART.type , 'CATE' )
                    p1 = addprop( obj , 'W_tr' );
                    p1.SetAccess = 'protected';
                    p2 = addprop( obj , 'ITT_tr' );
                    p2.SetAccess = 'protected';
                    if obj_optCART.honest
                        p3 = addprop( obj , 'W_est' );
                        p3.SetAccess = 'protected';
                        p4 = addprop( obj , 'ITT_est' );
                        p4.SetAccess = 'protected';
                        obj.W_tr  = ITT( ind_tr , : );
                        obj.W_est = ITT( ind_est , : );
                        obj.ITT_tr  = ITT( ind_tr , : );
                        obj.ITT_est = ITT( ind_est , : );
                    else
                        obj.W_tr    = ITT;
                        obj.ITT_tr  = ITT;
                    end
                
                % In case of fuzzy RDD only reset the value of ITT
                elseif strcmp( obj_optCART.type , 'CLATE' )
                    if obj_optCART.honest
                        obj.ITT_tr  = ITT( ind_tr , : );
                        obj.ITT_est = ITT( ind_est , : );
                    else
                        obj.ITT_tr  = ITT;
                    end
                end
                
                % If scaling was set for the running variable
                if obj_optCART.scale_runningvariable_w_c_2_zero
                    obj.X_tr = obj.X_tr - c;
                    if obj_optCART.honest
                        obj.X_est = obj.X_est - c;
                    end
                    obj.c = 0;
                end
            end
            % If scaling of the variables 
            if obj_optCART.scaleFeatures
                error('samples:InconsistentDataType','Scaling is not implemented at this moment...')
            end
            
            % If clustered SE is used
            if strcmp( obj_optCART.model , 'linear' )
                clV = iP.Results.cluster;
                if ~isempty( clV )
                    if all( [ n , 1 ] == size( clV ) ) && ( iscellstr( clV ) || isnumeric( clV ) )
                        if iscellstr( clV )
                            % If cellstring is given convert to numeric values
                            % for easier handling
                            [ aux , ~ , auxID ] = unique( clV );
                            cl_id0 = 1 : numel( aux );
                            clV  = cl_id0( auxID )';
                        end
                    else
                        error('samples:WrongInput:cluster' , ...
                            'Invalid type of data given for clusters, use numeric or cellstr vector, with the same number of observations and with one column!' )
                    end
                    if all( clV( 1 ) == clV )
                        error('samples:WrongInput:cluster' , ...
                            'Invalid input there is no variation in cluster, proper variable for clusters!' )
                    end
                    if any( isnan( clV ) )
                        error('samples:WrongInput:cluster' , ...
                            'There are NaN values in clusters! First clean your data!' )
                    end
                    if obj_optCART.honest
                        p1 = addprop( obj , 'cl_tr' );
                        p2 = addprop( obj , 'cl_est' );
                        p1.SetAccess = 'protected';
                        p2.SetAccess = 'protected';                        
                        obj.cl_tr  = clV( ind_tr , : );
                        obj.cl_est = clV( ind_est , : );
                    else
                        p1 = addprop( obj , 'cl_tr' );
                        p1.SetAccess = 'protected';
                        obj.cl_tr  = clV;
                    end
                end
            end
            update_sample( obj );
        end
        
        function val = get.n_tr( obj )
            val = size( obj.Y_tr , 1 );
        end
        function val = get.n_est( obj )
            if isprop( obj , 'Y_est' )
                val = size( obj.Y_est , 1 );
            else
                error('samples:InconsistentDataType','No such property as Y_est for sample object!')
            end
        end
        function val = get.K( obj )
            val = size( obj.Z_tr , 2 );
        end
        function val = get.L( obj )
            if isprop( obj , 'X_tr' )
                val = size( obj.X_tr , 2 );
            else
                error('samples:InconsistentDataType',...
                    'No covariates (X_tr) used during modelling the conditional expectations directly')
            end
        end
    end
    
    %% Create Cross-Validation samples
    methods
        function update_sample( obj )
            if isprop( obj , 'X_tr' )
                obj.log_X_tr   = true;
            end
            if isprop( obj , 'W_tr' )
                obj.log_W_tr   = true;
            end
            
            if isprop( obj , 'ITT_tr' )
                obj.log_ITT_tr = true;
            end
            
            if isprop( obj , 'cl_tr' )
                obj.log_cl_tr  = true;
            end
            
        end
        [ train_sample , test_sample ] = create_CV_samples( obj_samples , obj_optCART , K );
    end
    
    
end