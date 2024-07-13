function [X1 , X2] = solver_vhl(y, Bs, s, n)
 
% K = 2


if mod(n,2) == 0
    n1 = n/2;
else
    n1 = (n+1)/2;
end
n2 = n+1 - n1;

cvx_begin quiet
    %cvx_precision best
    variable X1(s,n) complex;
    variable X2(s,n) complex;
    
    expression obj1(s*n1, n2);
    expression obj2(s*n1, n2);
    % hankel
    for j1 = 1:n1
        for j2 = 1:n2
            row_idx1 = (j1-1)*s+1:j1*s;
            obj1(row_idx1, j2) = X1(:, j1+j2-1);
            row_idx2 = (j1-1)*s+1:j1*s;
            obj2(row_idx2, j2) = X2(:, j1+j2-1);
            
        end
    end
    
    minimize norm_nuc(obj1)+norm_nuc(obj2);
    subject to
        % measurement
        
        y == diag(Bs(:,:,1)*X1) + diag(Bs(:,:,2)*X2);
        
    
    
cvx_end 




end
