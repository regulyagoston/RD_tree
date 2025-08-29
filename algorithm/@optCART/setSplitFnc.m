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
        % Method implemented with parametric RDD
        obj.split_fnc = @split_CT_RDD;
    case 'causal-tree-npRDD-uni'
        % Method implemented for RDD with unified bandwidth
        obj.split_fnc = @split_CT_npRDD_uni;
    case 'causal-tree-npRDD-lbl'
        % Method implemented for RDD with leaf-by-leaf bandwidths
        obj.split_fnc = @split_CT_npRDD_lbl;
    otherwise
        error( 'No such splitting type implemented!' )
end

end