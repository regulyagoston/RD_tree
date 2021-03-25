%% Create polinomials for X

function XX = createPolinomials( X , orderPoly )

[ N , M ] = size( X );

if M ~= 1
   error('X must be a Nx1 vector and not a matrix!') 
end

if orderPoly == 1
    XX = X;
    return
else
    XX = NaN( N , orderPoly );
    for i = 1 : orderPoly
        XX( : , i ) = X.^i;
    end
end

end