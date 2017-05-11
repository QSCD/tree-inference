function densityAfterDivision =  applyCellDivision_matrix(density,statespace)

assert(length(statespace)+1==length(density))

maxState=floor(statespace(end)./2);
max_consideredStates = 0:maxState;

matrix = zeros(length(statespace));
for i =1:length(max_consideredStates)
    
    currentState =  max_consideredStates(i);
    weights = zeros(1,length(statespace));
    if currentState == 0 % only has input from above, not below
        ix1 = 2.*currentState ==statespace;
        ix2 = 2.*currentState+1 ==statespace;
        
        weights(ix1) =1;
        weights(ix2) = 0.5;
        
        
    elseif maxState==currentState
        ix1 = 2.*currentState ==statespace;
        ix2 = 2.*currentState-1 ==statespace;
        
        weights(ix1) =1;
        weights(ix2) = 0.5;
        
    else
        ix1 = 2.*currentState ==statespace;
        ix2 = 2.*currentState+1 ==statespace;
        ix3 = 2.*currentState-1 ==statespace;
        weights(ix1) =1;
        weights(ix2) = 0.5;
        weights(ix3) = 0.5;
    end
    matrix(i,:)=weights;
    
end
%add the rows of the absorbing state
%it does not gain anything and does not conribute to anything hwnce just 0
matrix = [matrix zeros(size(matrix,1),1)];
matrix = [matrix; zeros(1,size(matrix,2))];
matrix = sparse(matrix);
densityAfterDivision = density*matrix';  % more standard : (matrix * dens')'
end