classdef paramLP < handle
    
    %% Order of Polinomial and its derivative
    properties
        p      = 1;
        deriv  = 0;
    end
    
    %% Estimation options
    properties
        eval   = [];
        neval  = 30;
        subset = [];
    end
    
    %% Kernel Properties
    properties
        h
        kernel   = 'Triang';
        bwselect = 'mse-dpi';
        bwcheck  = 21;
        bwregul  = true;
        imsegrid = 30;
        interior = false;
        all      = false;       %If all bandwidth is estimated
    end
    
    %% Variance Properties
    properties
        vce     = 'hc0';
        cluster = [];
        nnmatch = 3;
        level   = 0.95;
    end
    
    %% Bias Correction Properties
    properties
        q   = 2;
        b
        rho = 1;
    end
    
    %% RDD Properties
    properties
       scalepar   = 1;      % Scaling factor for RD parameter of interest 
       scaleregul = 1;      % Scaling factor for the regularization for BW selection
       sharpbw    = false;  % Perform fuzzy RD estimation using bw for sharp RD model
    end
    
    methods
        function obj = paramLP()
            return;
        end
        
        function set.p( obj , val )
            if ~( numel( val ) == 1 && mod( val , 1 ) == 0 && val > 0 )
                error( 'paramLP:WrongInput:p' ,'Invalid value for order of local polinomial, use a positive integer!')
            elseif ~( isempty( obj.deriv ) || val > obj.deriv )
                error( 'paramLP:WrongInput:p' ,'Invalid value for order of local polinomial, use higher value than the reqiured derivative!')
            else
                obj.p = val;
            end
        end
        
        function set.deriv( obj , val )
            if ~( numel( val ) == 1 && mod( val , 1 ) == 0 && val >= 0 )
                error( 'paramLP:WrongInput:deriv' ,'Invalid value for the derivative, use zero or a positive integer!')
            else
                if obj.p <= val
                    warning( 'Order of polinomial has been set to higher value due to higher value of derivative.');
                    obj.p = val;
                end
                if obj.q <= obj. p
                    obj.q = obj.p + 1;
                end
                obj.deriv = val;
            end
        end
        
        function set.q( obj , val )
            if ~( numel( val ) == 1 && mod( val , 1 ) == 0 && val > 0 )
                error( 'paramLP:WrongInput:q' ,'Invalid value for the auxillary order of local polinomial, use a positive integer!')
            else
                obj.q = val;
            end
        end
        
        function set.nnmatch( obj , val )
            if ~( numel( val ) == 1 && mod( val , 1 ) == 0 && val > 0 )
                error( 'paramLP:WrongInput:nnmatch' ,'Invalid value for number of neighbours to match, use a positive integer!')
            else
                obj.nnmatch = val;
            end
        end
        
        function set.bwselect( obj , val )
           if ~ischar( val ) && ~any( strcmpi( val , {'imse-dpi','mse-dpi','mse-rd','manual'} ) )
               error( 'paramLP:WrongInput:bwselect' , ...
                   ['Invalid value for bandwidth selector, use one of the following: \n',...
                    'imse-dpi: integrateg mean squared error based on direct plug in \n' ,...
                    'mse-dpi : mean squared error based on direct plug in \n',...
                    'mse-rd  : mean squared error for RD \n' , ...
                    'range   : use all the observation available \n' , ...
                    'manual  : manually given value for bandwidth'])
           else
               obj.bwselect = val;
           end
        end
        
        function set.vce( obj , val )
            if ~ischar( val ) && ~any( strcmpi( val , {'nn','hc0','hc1','hc2','hc3','cluster','nncluster'} ) )
               error( 'paramLP:WrongInput:vce' , ...
                   ['Invalid value for variance-covariance estimator, use one of the following: \n',...
                    'nn: nearest neighbour estimator \n',...
                    'hc0-hc3: different versions of heteroscedastic estimator \n',...
                    'cluster\nncluster: clustered varaiance' ] )
           else
               obj.vce = val;
            end
        end
        
        function set.kernel( obj , val )
            if ~any( strcmp( val , { 'Uni' , 'Epa' , 'Norm' , 'Triang' } ) )
                error( 'paramLP:WrongInput:kernel' , ...
                    ['Invalid input for kernel density, use one of the following: \n',...
                    'Uni    - uniform \n' , ...
                    'Epa    - epanechnikov \n' , ...
                    'Norm   - normal \n' , ...
                    'Triang - triangular kernel '] )
            else
                obj.kernel = val;
            end
        end
        
        function set.eval( obj , val )
            if ~( isempty( val ) || ( isnumeric( val ) && ~any( isnan( val ) ) ) )
                error( 'paramLP:WrongInput:eval' , 'Invalid input for evaluation vector, use numeric, Kx1 vector!')
            else
                obj.eval = val( : );
            end
        end
        function set.neval( obj , val )
            if ~( isempty( val ) || ( numel( val ) == 1 && mod( val , 1 ) == 0 && val > 0 ) )
                error( 'paramLP:WrongInpit:neval' , 'Invalid input for number of evaluation points, use positive integer!')
            else
                obj.neval = val;
            end
        end
        function set.level( obj , val )
            if ~( isnumeric( val ) && numel( val ) == 1 && val > 0 && val < 1 )
                error( 'paramLP:WrongInput:level' , 'Invalid value for confidence levels, use a numeric value between 0 and 1')
            else
                obj.level = val;
            end
        end
%         function set.rho( obj , val )
%             if ~( numel( val ) == 1 && val > 1 )
%                 error( 'paramLP:WrongInput:rho' , 'Invalid value for rho, use a value greater than 1!' )
%             else
%                 obj.rho = val;
%             end
%         end
    end
    
    methods
        function clearEstimates( obj )
            obj.eval = [];
            % obj.neval = 30;
            obj.subset = [];
            obj.h = [];
            obj.cluster = [];
            obj.b = [];
        end
        
    end
    
end