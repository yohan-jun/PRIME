function U = svd_left(A, r)
% Left singular matrix of SVD (U matrix)
% parameters: matrix, rank (optional)
if nargin < 2
    [U,E] = eig(A*A'); % surprisingly, it turns out that this is generally faster than MATLAB's svd, svds, or eigs commands
    [~,idx] = sort(abs(diag(E)),'descend');
    U = U(:,idx);
else
    [U,~] = eigs(double(A*A'), r);
    U = cast(U, 'like', A);
end
end
