function  gstz = Gstar(Z,D, s)
[N,n2] = size(Z);
n1 = N/s;

gstz = zeros(s,n1+n2-1);
for i = 1:(n1+n2-1)
    if i <= n1
        comp = zeros(s,1);
        for j = 1:i
            comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
        end
        gstz(:,i) = 1/D(i)*comp;
    else
        comp = zeros(s,1);
        for j = i+1-n1:n2
            comp = comp + Z(s*(i-j)+1:s*(i+1-j),j);
        end
        gstz(:,i) = 1/D(i)*comp;
    end
end
end