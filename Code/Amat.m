function matA = Amat(z, L)
    A = zeros(L+1);
    A(1,1) = 1;

    for q = 2:L+1
        A(1,q) = 0;
    end

    for q = 1:L
        A(2,q) = z(q);
    end

    for p = 3:L+1
        for q = 1:L+1-p+1
            A(p,q) = A(p-1,1)*A(p-2,q+1) - A(p-2,1)*A(p-1,q+1);
        end
    end
    
    matA = num2cell(A);
    
end