function isposdef = ismatrix_positive(A)
% 判断矩阵是否正定 如果 p = 0，则输入矩阵是对称正定矩阵
    % Cholesky Decomposition
    [~, p] = chol(A);
    if p == 0
        isposdef = true;  % A is positive definite
    else
        isposdef = false; % A is not positive definite
    end
end