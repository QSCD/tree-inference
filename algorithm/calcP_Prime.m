% calculates the "commonTree equaivalent" (that is, the undifferentiated part of the tree) in the speedy version of the likelihood
% just considering the branch above the differentiation point and correcting for the double counting
% of cells by only using left  sisters
%
% returns PPrime, where PPrime(i) corresponds to the undiff prob (primed) of cell(i) 
% in other words, cells is just the corresponding ordering of PPrime

% e.g if we want the the prob of a tree being undifferentiated at cells 2,6,7 (make sure this is an entire set of leaves no open ends)
% its just PPrime(cell2)* PPrime(cell6) * PPrime(cell7): no need to worry about counting things
% twice, its already accounted for
%
%
% relevantCells_perEntry{i}: the cells that have been used to calc PPrime(i)
%           if thrown together, it should give no overlap!!
%
% p_undiff_perCell: the raw accumulated hazard for this cell not to diff
function [PPrime,cells,p_undiff_perCell] = calcP_Prime(precalcStruct,k_diff,currentTransformed)

diff_offset = k_diff(1);
diff_slope = k_diff(2);

% accounting for the unoberserved stuff before the movie
% acutally this should not be necessary, because we dont have to cut the trees any more
% ===================
% NOTICE: THIS IS MULTIPLIED ONTO THE WHOLE STUFF WHEN WE SELECT THE ROOT OF THE TREE
if isfield(currentTransformed, 'startingTime')
    movieStart = currentTransformed.commonUpTreeCorrection_start_stop(1); % notice that this is not necces.  movie start (-> massive overlap if the trees com from the same huge tree), but the time until it merges with some other subtree
    %note that the .commonUpTreeCorrection_start_stop(1) should be pretty close to
    %startingTime!!
    startingTime = currentTransformed.startingTime;
    additionalFactor = exp(-(diff_offset.*(startingTime-movieStart)+ 0.5.*diff_slope .* (startingTime.^2-movieStart.^2)));

    if additionalFactor~=1 || startingTime~=0
        warning('should not do that!!: dsont use the faster likelihood with the shifted starting times ')
    end
else
    additionalFactor=1;
    startingTime = 0;
end
    
    
%lets calc the hazard for each cel lfirst, 
%nothing special here!
cells = precalcStruct.cells;

currentBirth = precalcStruct.birthtimes + startingTime;
currentDeath = precalcStruct.deathtimes + startingTime;
p_undiff_perCell= exp(-(diff_offset.*(currentDeath-currentBirth)+ 0.5.*diff_slope .* (currentDeath.^2-currentBirth.^2)));


% now that we have it for each cell, combinethem into PPrime(i)
% which is the prob to be undiff until (and including )i
% however accounting only for left branches
% that is, calc the prb on the rootpath until the current cell is the right sister somwhere
% (however, we also have to check the the left sister indeed exists!)
PPrime = zeros(1,length(cells));
for i = 1:length(cells)

   cIx = precalcStruct.relevantCells_perEntry_IX{i}; % we precalculated which entries of p_undiff_perCell are used to calc PPrime(i)!
   if ~precalcStruct.isRootFlag(i)  %it the current cells entry does not stretch to the root
       PPrime(i) = prod(p_undiff_perCell(cIx)); % equals 1 if no cell is relevant (meaning were at the root), but thats just fine!
   else
       PPrime(i) = prod(p_undiff_perCell(cIx)).*additionalFactor; %for th root cell we also include the additional factor the accounts for the "unobserved" cells before the root
   end
end