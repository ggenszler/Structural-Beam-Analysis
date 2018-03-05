% function solves degree of static indeterminacy 
function degindet=indet(y)
    rxns=0;
    [m,n]=size(y);
    for x=1:m
        rxns=rxns+3-y(x,1);
    end
    % for every row in the matrix y, the first column's entry is examined
    % and used to calculate the number of reactions
    degindet=rxns-3;
    % degree indeterminate is equal to the number of reactions minus thrice
    % the number of sections. However, the number of sections will always
    % be one for the cases this code can handle.
end