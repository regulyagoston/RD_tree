%% Set the splitting function for CART
function setSplitFnc( obj )

switch obj.splitType
    case 'pctg_wM'
        % percentile with median enforced
        obj.split_fnc = @split_pctg_wM;
    case 'pctg'
        % percentile without median enforced
        obj.split_fnc = @( Zi , nS )( prctile( Zi , linspace( 0 , 100 , nS ) ) );
    case 'uni_support'
        % Equally distanced
        obj.split_fnc = @( Zi , nS )( linspace( min( Zi ) , max( Zi ) , nS ) );
    case 'rpart'
        % Method implemented in rpart
        obj.split_fnc = @( Zi , a )( unique( Zi ) );
    case 'causal-tree-AI'
        % Method implemented in causal-trees
        obj.split_fnc = @splitAthey_Imbens;
    case 'causal-tree-RDD'
        % Method implemented in causal-trees
        obj.split_fnc = @splitRDD;
    otherwise
        error( 'No such splitting type implemented!' )
end

end