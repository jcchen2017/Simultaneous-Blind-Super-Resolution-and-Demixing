function GX = G(X)

[s,n]  = size(X);
if mod(n,2) == 0
    n1 = n/2;
    DD = [1:n1 n1 n1-1:-1:1].';
else
    n1 = (n+1)/2;
    DD = [1:n1 n1-1:-1:1].';
end
n2 = n+1 - n1;
D = sqrt(DD);

X_scaled = X*diag(1./D);
GX = zeros(s*n1, n2);
for j1 = 1:n1
    for j2 = 1:n2
        row_idx = (j1-1)*s+1:j1*s;
        GX(row_idx, j2) = X_scaled(:, j1+j2-1);
    end
end
end