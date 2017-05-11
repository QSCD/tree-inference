%essentially, extracts the S from the D defined via the recursion
% D(i) = S(i) + D(2i)D(2i+1)
function marginals = D2marginal(D_array,cellNrs)

assert(length(D_array)==length(cellNrs))

marginals = zeros(length(cellNrs),1);
for i =1:length(cellNrs)
    cC = cellNrs(i);
    lD = cC.*2;
    rD = lD+1;
    
    ix_lD = cellNrs==lD;
    ix_rD = cellNrs==rD;
    
    if all(ix_lD==0)
       ld_contrib = 1; 
    else
        ld_contrib = D_array(ix_lD);
    end
    
    if all(ix_rD==0)
       rd_contrib = 1; 
    else
        rd_contrib = D_array(ix_rD);
    end
    
    if all(ix_rD==0) & all(ix_lD==0) % it th4e current cell is a leave, we dont have to subtract contribs from daughters
         marginals(i) = D_array(i);
    else %otherwise we have to subtract their contribs
        marginals(i) = D_array(i) - ld_contrib*rd_contrib;
    end
    
end