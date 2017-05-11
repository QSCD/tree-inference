% this matrix M , applied to a density row vector x simulates the cell
% division:
% divided = x*M 
% note that the matrix has to be applied from the right!!
function divMatrix = calcDivisionMatrix(statespace)


divMatrix = zeros(length(statespace));

PERFECT_SYMMETRY =1;
if PERFECT_SYMMETRY == 1
    maxState=ceil(statespace(end)./2); %notre that this is diferent when doing perfectly symetric or not: if symmetric we have to round up
    max_consideredStates = 0:maxState;    
    for i =1:length(max_consideredStates)
        currentState =  max_consideredStates(i);
        weights = zeros(1,length(statespace));
        if currentState == 0 % 0 has only input from 0 (state one divides into (1,1) already)
            ix1 = currentState ==statespace;
            weights(ix1) =1;
            
        elseif maxState==currentState
            if mod(currentState,2) == 0 %if its even only the double state contributres
                ix1 = 2.*currentState ==statespace;               
                ix2 = 2.*currentState-1 ==statespace;
                weights(ix1) =1;
                weights(ix2) =1;
            else % the double and the double +1 contribute
                ix1 = 2.*currentState ==statespace;
                ix2 = 2.*currentState-1 ==statespace;
                
                weights(ix1) =1;
                weights(ix2) =1;
            end
            
        else
            ix1 = 2.*currentState ==statespace;
            ix2 = 2.*currentState-1 ==statespace;
            
            weights(ix1) =1;
            weights(ix2) = 1;
        end
        divMatrix(i,:)=weights;
    end
else
    maxState=floor(statespace(end)./2);
    max_consideredStates = 0:maxState;
    for i =1:length(max_consideredStates)

        currentState =  max_consideredStates(i);
        weights = zeros(1,length(statespace));
        if currentState == 0 % only has input from above, not below
            ix1 = 2.*currentState ==statespace;
            ix2 = 2.*currentState+1 ==statespace;

            weights(ix1) =1;
            weights(ix2) = 0.5;


        elseif maxState==currentState %wenn man am anderen ende des statespace ist, kommt kein betrag von größeren states
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
        divMatrix(i,:)=weights;

    end
end
assert(all(sum(divMatrix,1)==1),'make sure that it sums to one')
%add the rows of the absorbing state
%it does not gain anything and does not conribute to anything hwnce just 0
divMatrix = [divMatrix zeros(size(divMatrix,1),1)];
divMatrix = [divMatrix; zeros(1,size(divMatrix,2))];
divMatrix = sparse(divMatrix');

end