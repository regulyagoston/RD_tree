


function [ pred , pred_var ] = predict_forest( forest , inbag , Z )

nF = numel( forest );
N  = size( Z , 1 ); 
pred_ij = NaN( N , nF );

for i = 1 : nF
    pred_ij( : , i ) = predict( forest( i ) , Z );
end

% Make prediction from pred_ij and if needed get the variances
if nargout == 1
    pred = mean( pred_ij , 2 );
else
    [ pred , pred_var ] = forestInfJack( pred_ij , inbag , true , 1 : nF );
end


end