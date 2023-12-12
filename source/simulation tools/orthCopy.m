function Q = orthCopy(X)
    % Modified Gram-Schmidt orthonormalization 
    % which produces the same result as [Q,R]=qr(X,0)
    % Written by Mo Chen (sth4nth@gmail.com).
    [d,n] = size(X);
    m = min(d,n);
    R = zeros(m,n);
    Q = zeros(d,m);
    for i = 1:m
        v = X(:,i);
        for j = 1:i-1
            R(j,i) = Q(:,j)'*v;
            v = v-R(j,i)*Q(:,j);
        end
        R(i,i) = norm(v);
        Q(:,i) = v/R(i,i);
    end
end