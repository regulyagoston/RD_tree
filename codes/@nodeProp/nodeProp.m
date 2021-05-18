%% This value class contains the essential informations for a given node in the tree
classdef nodeProp 
    %% Informations about the node
    properties
        Checked  = false;   % Whether node has been checked for further splits
        Failed   = 0;       % Failed variables in the nodes, and why:
                            %   0 - not failed
                            %   1 - 'No Increase in measurement criteria'
                            %   2 - 'No splitting value remained'
                            %   3 - 'Reached maximum depth (level) for the tree'
                            %   4 - 'Reached maximum number of nodes'
                            %   5 - 'Reached maximum number of leaves'
    end
    %% CORE PROPERTIES: Subset selection
    properties
        logID_tr            % Logicals which select the sample for the leaf - TRAINING SAMPLE
        subsetIdx           % Which features from Z has been used (index)
        splitVal            % Where is the split value from Z( : , subsetIdx )
        LEq                 % Less or Equal  - contains the less or equal part of the sample (true)
        feature_chck        % Check whether feature is a candidate split
                            %   true  - it searches along the feature
                            %   false - it neglects the feature from this
                            %           node down the branch
    end
    %% Estimated values - structure
    properties
       est_struct           % Estimated values in a structure
    end
    %% Simple (adaptive) estimation
    properties
        n_j_tr      % Number of observations for the TRAINING SAMPLE
    end
    %% Honest estimation
    properties
        logID_est   % Logicals which select the sample which is member of the leaf - ESTIMATION
        n_j_est     % Number of observations for the ESTIMATION SAMPLE
    end
    
    methods
        function obj = nodeProp( subsetIdx , splitVal , LEq , Failed , feature_chck )
            if nargin == 0
                return;
            end
            % Node with all data
            obj.subsetIdx    = subsetIdx;
            obj.splitVal     = splitVal;
            obj.LEq          = LEq;
            obj.Failed       = Failed;
            obj.feature_chck = feature_chck;
        end
    end
    
    methods
        % Update the Node with estimated values
        function obj = updateNode( obj , data , optTree , tr )
            % Update basic values
            if isempty( obj.logID_tr )
                obj.logID_tr  = get_logIdx( obj , data.Z_tr  );
            end
            obj.n_j_tr = sum( obj.logID_tr );
            if isempty( obj.logID_est ) && optTree.honest
                obj.logID_est  = get_logIdx( obj , data.Z_est  );
                obj.n_j_est = sum( obj.logID_est );
            end
            
            %% If Leaf-by-leaf criterion -> estimate the relevant values
            if strcmp( optTree.eval_level , 'leaf-by-leaf' )
                msg1 = false;
                % Treatment estimations for CATE
                if strcmp( optTree.type , 'CATE' )
                    if strcmp( optTree.model , 'mean' )
                        obj.est_struct = est_treat_lbl_mean( obj , data , optTree , tr );
                    elseif strcmp( optTree.model , 'linear' )
                        obj.est_struct = est_treat_lbl_OLS( obj , data , optTree , tr );
                    else
                        msg1 = true;
                    end
                % Treatment estimations for CLATE
                elseif strcmp( optTree.type , 'CLATE' )
                    if strcmp( optTree.model , 'linear' )
                        obj.est_struct = est_treat_lbl_OLS( obj , data , optTree , tr );
                    else
                        msg1 = true;
                    end
                % Prediction estimation
                elseif strcmp( optTree.type , 'prediction' )
                     if strcmp( optTree.model , 'mean' )
                         % In case of simple mean, just get the data and
                         %     estimate the means and replace observations with
                         %     the mean
                         Y                      = get_node_tr_Sample( obj , data );
                         mY                     = mean( Y );
                         obj.est_struct.Y_hat_i = repmat( mY , [ obj.n_j_tr , 1 ] );
                     else
                         msg1 = true;
                     end
                else
                    msg1 = true;
                end
                if msg1 
                    error('nodeProp:updateNode:incorrectType',...
                           'No such updating function defined, while updating the node!')
                end
            end
        end
        % Predict a given sample with an already estimated node
        % Update the Node with estimated values
        function obj = predictNode( obj , obj_tr , data , optTree , tr )
            % Update basic values
            if isempty( obj.logID_tr )
                obj.logID_tr  = get_logIdx( obj , data.Z_tr  );
            end
            obj.n_j_tr = sum( obj.logID_tr );
            if isempty( obj.logID_est ) && optTree.honest
                obj.logID_est  = get_logIdx( obj , data.Z_est  );
                obj.n_j_est = sum( obj.logID_est );
            end
            
            %% If Leaf-by-leaf criterion -> estimate the relevant values
            if strcmp( optTree.eval_level , 'leaf-by-leaf' )
                msg1 = false;
                % CATE estimations
                if strcmp( optTree.type , 'CATE' )
                    if any( strcmp( optTree.model , {'mean','linear'} ) )
                        obj.est_struct = obj_tr.est_struct;
                        tau_i = obj_tr.est_struct.tau_i( 1 );
                        if tr
                            obj.est_struct.tau_i = repmat( tau_i , [ obj.n_j_tr , 1 ] );
                        else
                            if obj.n_j_est == 0 
                                obj.est_struct.tau_i = tau_i;
                            else
                                obj.est_struct.tau_i = repmat( tau_i , [ obj.n_j_est , 1 ] );
                            end
                        end
                    else
                        msg1 = true;
                    end
                % CATE estimations
                elseif strcmp( optTree.type , 'CLATE' )
                    if strcmp( optTree.model , 'linear' )
                        obj.est_struct = obj_tr.est_struct;
                        tau_i = obj_tr.est_struct.tau_i( 1 );
                        tau_y_i = obj_tr.est_struct.tau_y_i( 1 );
                        tau_w_i = obj_tr.est_struct.tau_w_i( 1 );
                        if tr
                            obj.est_struct.tau_i   = repmat( tau_i , [ obj.n_j_tr , 1 ] );
                            obj.est_struct.tau_y_i = repmat( tau_y_i , [ obj.n_j_tr , 1 ] );
                            obj.est_struct.tau_w_i = repmat( tau_w_i , [ obj.n_j_tr , 1 ] );
                        else
                            if obj.n_j_est == 0 
                                obj.est_struct.tau_i   = tau_i;
                                obj.est_struct.tau_y_i = tau_y_i;
                                obj.est_struct.tau_w_i = tau_w_i;
                            else
                                obj.est_struct.tau_i   = repmat( tau_i , [ obj.n_j_est , 1 ] );
                                obj.est_struct.tau_y_i = repmat( tau_y_i , [ obj.n_j_est , 1 ] );
                                obj.est_struct.tau_w_i = repmat( tau_w_i , [ obj.n_j_est , 1 ] );
                            end
                        end
                    else
                        msg1 = true;
                    end
                % Prediction estimation
                elseif strcmp( optTree.type , 'prediction' )
                     if strcmp( optTree.model , 'mean' )
                         % In case of simple mean, just get the data and
                         %     estimate the means and replace observations with
                         %     the mean
                         mY = obj_tr.est_struct.Y_hat_i( 1 );
                         obj.est_struct.Y_hat_i = repmat( mY , [ obj.n_j_tr , 1 ] );
                     else
                         msg1 = true;
                     end
                else
                    msg1 = true;
                end
                if msg1 
                    error('nodeProp:updateNode:incorrectType',...
                           'No such updating function defined, while updating the node!')
                end
            end
        end
        % Logical Indexing for training sample
        function logId = get_logIdx( obj , Z )
            n = size( Z , 1 );
            if ~isempty( obj.subsetIdx )
                nC = numel( obj.subsetIdx );
                logId_a = false( [ n , nC ] );
                for i = 1 : nC
                    if obj.LEq( i )
                        logId_a( Z( : , obj.subsetIdx( i ) ) <= obj.splitVal( i ) , i ) = true;
                    else
                        logId_a( Z( : , obj.subsetIdx( i ) ) >  obj.splitVal( i ) , i ) = true;
                    end
                end
                logId = all( logId_a , 2 );
            else
                logId = true( [ n , 1 ] );
            end
        end
        % Get training sample
        function [ Y , X , W , ITT , cl ] = get_node_tr_Sample( obj , data )
            if isempty( obj.logID_tr )
                obj.logID_tr  = get_logIdx( obj , data.Z_tr  );
            end
            Y = data.Y_tr( obj.logID_tr , : );
            if data.log_X_tr
                X = data.X_tr( obj.logID_tr , : );
            else
                X = NaN;
            end
            if data.log_W_tr
                W = data.W_tr( obj.logID_tr , : );
            else
                W = NaN;
            end
            if data.log_ITT_tr
                ITT = data.ITT_tr( obj.logID_tr , : );
            else
                ITT = NaN;
            end
            if data.log_cl_tr
                cl = data.cl_tr( obj.logID_tr , : );
            else
                cl = NaN;
            end
        end
        % Get estimation sample
        function [ Y , X , W , ITT , cl ] = get_node_est_Sample( obj , data )
            if isempty( obj.logID_est )
                obj.logID_est  = get_logIdx( obj , data.Z_est  );
            end
            Y = data.Y_est( obj.logID_est , : );
            if isprop( data , 'X_est' )
                X = data.X_est( obj.logID_est , : );
            else
                X = NaN;
            end
            if isprop( data , 'W_est' )
                W = data.W_est( obj.logID_est , : );
            else
                W = NaN;
            end
            if isprop( data , 'ITT_est' )
                ITT = data.ITT_est( obj.logID_est , : );
            else
                ITT = NaN;
            end
            if isprop( data , 'cl_est' )
                cl = data.cl_est( obj.logID_est , : );
            else
                cl = NaN;
            end
        end
        % Get features which is not failed
        function [ idx , n_id ] = get_valid_features( obj )
            idx  = find( obj.feature_chck );
            n_id = numel( idx );
        end
    end
    
    methods
        %% Clear all the estimated values from the node (only core properties are inherited)
        function obj2 = clearNode( obj )
            % Copy only the split/tree information about the node
            obj2 = nodeProp( obj.subsetIdx , obj.splitVal , obj.LEq , obj.Failed , obj.feature_chck );
        end
        function obj2 = reset_NodeIdx( obj , data )
            % Copy only the split/tree information about the node
            obj2 = nodeProp( obj.subsetIdx , obj.splitVal , obj.LEq , obj.Failed , obj.feature_chck );
            obj2.est_struct = obj.est_struct;
            obj2.logID_tr   = get_logIdx( obj2 , data.Z_tr  );
            obj2.n_j_tr     = sum( obj2.logID_tr );
            if optTree.honest
                obj2.logID_est = get_logIdx( obj2 , data.Z_est  );
                obj2.n_j_est   = sum( obj2.logID_est );
            end
        end
        %% Check whether the two nodes are the same between the tolerance range (only core properties)
        function tf = issameNode( node1 , node2 , tolerance )
            tf = false;
            if all( node1.subsetIdx == node2.subsetIdx ) && ...
               all( node1.LEq == node2.LEq ) && ...
               all( node1.Failed == node2.Failed )
                % Check for the split values
                nS = numel( node1.subsetIdx );
                chckSplit = false( nS , 1 );
                for s = 1 : nS 
                    logC = node1.splitVal( s ) - tolerance <= node2.splitVal( s ) && ...
                           node1.splitVal( s ) + tolerance >= node2.splitVal( s );
                    chckSplit( s ) = logC;
                end
                if all( chckSplit )
                    tf = true;
                end
            end
        end
    end
end

