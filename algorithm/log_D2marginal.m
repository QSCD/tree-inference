%essentially, extracts the S from the D defined via the recursion
% D(i) = S(i) + D(2i)D(2i+1)
%
% this works in logspace
function log_marginals = log_D2marginal(logD_array,cellNrs)

assert(length(logD_array)==length(cellNrs))

log_marginals = zeros(length(cellNrs),1);
for i =1:length(cellNrs)
    cC = cellNrs(i);
    lD = cC.*2;
    rD = lD+1;
    
    ix_lD = cellNrs==lD;
    ix_rD = cellNrs==rD;
    
    if all(ix_lD==0)
       log_ld_contrib = 0;  % ld_contrib = 1; <-lin space
    else
        log_ld_contrib = logD_array(ix_lD);
    end
    
    if all(ix_rD==0)
       log_rd_contrib = 0; %rd_contrib = 1; 
    else
        log_rd_contrib = logD_array(ix_rD);
    end
    
    if all(ix_rD==0) & all(ix_lD==0) % it th4e current cell is a leave, we dont have to subtract contribs from daughters
         log_marginals(i) = logD_array(i);
    else %otherwise we have to subtract their contribs (in linspace)
        % i.e marginal = D(i)-leftContr*rightContr
        %  logmarginal = log[ exp(logD(i)) - exp(logleft+logright)]
        log_marginals(i) = LogDiffExp(logD_array(i),  log_rd_contrib+log_ld_contrib);
    end
    
end
end



% logExpSum trick but for two varialbes and difference instead of sum
function overall = LogDiffExp(X1,X2)
maxX = max(X1,X2);

overall = maxX + log( exp(X1-maxX) - exp(X2-maxX) );
end