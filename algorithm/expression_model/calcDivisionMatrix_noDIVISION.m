% this matrix M , applied to a density row vector x simulates the cell
% division:
% divided = x*M 
% note that the matrix has to be applied from the right!!

%NOTE THIS ONE HERE DOES ACTUAL NOT DIVIDE THE PROTEINS BUT KEEPS THEM
function divMatrix = calcDivisionMatrix_noDIVISION(statespace)


divMatrix = eye(length(statespace));

assert(all(sum(divMatrix,1)==1),'make sure that it sums to one')
%add the rows of the absorbing state
%it does not gain anything and does not conribute to anything hwnce just 0
divMatrix = [divMatrix zeros(size(divMatrix,1),1)];
divMatrix = [divMatrix; zeros(1,size(divMatrix,2))];
divMatrix = sparse(divMatrix');

end