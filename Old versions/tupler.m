function [tup, j, ell] = tupler(tup,j,ell,N)
%used to recursively find all tuples which sum to N

    %to stop infinite looping of infinite zeros:
    if ell==N
        tup{j}(N)=N-sum(tup{j});
        return;
    end
    
    %main recursive bit:
    for n=0:N
        ell=ell+1;
        tup{j}(ell)=n;
        if sum(tup{j})==N
           %start a new tuple
           j=j+1; 
           ell=0;
           return; 
        else
            [tup, j, ell] = tupler(tup,j,ell,N);
            if sum(tup{j})==N
                %start a new tuple
                j=j+1; 
                ell=0;
                return;
            end
        end
    end
end