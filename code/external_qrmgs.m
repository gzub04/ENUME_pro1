% it's the same function as the one in the class, but modified so it can
% accept some other matrices not only the ones already inside the class
function [Q, R] = external_qrmgs(D)
    %QR (thin) factorization using modified Gram-Schmidt algorithm
    %for full rank real-valued and complex-valued matrices
    [m, n] = size(D);
    d = zeros(1, n);
    Q = zeros(m, n);
    R = zeros(n, n);
    %factorization with orthogonal (not orthonormal) columns of Q:
    for i = 1 : n
        Q(:, i) = D(:, i);
        R(i, i) = 1;
        d(i) = Q(:, i)' * Q(:, i); 
        for j = (i + 1) : n
            R(i, j) = (Q(:, i)' * D(:, j)) / d(i);
            D(:, j) = D(:, j) - R(i, j) * Q(:, i);
        end
    end
    %column normalization (columns of Q orthonormal):
    for i = 1 : n
        dd = norm(Q(:, i));
        Q(:, i) = Q(:, i) / dd;
        R(i, i:n) = R(i, i:n) * dd;
    end
end

