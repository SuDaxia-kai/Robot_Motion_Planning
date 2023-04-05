function Q = getQ(n_seg, n_order, ts, r)
    Q = [];
    T = zeros((n_order-r)*2+1,1);
    for k = 1:n_seg
%         Q_k = [];
        %#####################################################
        % STEP 1.1: calculate Q_k of the k-th segment 
        for i = 1:(n_order-r)*2+1
            T(i) = ts(k)^i;
        end
        Q_k = zeros(n_order+1, n_order+1);
        for row = r+1:n_order+1
            for line = row:n_order+1
                k1 = row-r-1;
                k2 = line-r-1;
                k3 = k1+k2+1;
                Q_k(row,line) = prod(k1+1:k1+r)*prod(k2+1:k2+r)/k3*T(k3);
                % 由于Q是对称矩阵，所以只需要算出一半，另外一般对称过去就行
                Q_k(line,row) = Q_k(row,line);
            end
        end
        Q = blkdiag(Q, Q_k);
    end
end