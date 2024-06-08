function [T,Q] = lanczos(A,B,iterMax)

[n,nB] = size(B);
Q = zeros(n,n);
alpha = cell(n,1);
beta = cell(n,1);
B = orth(B);
for j = 1:iterMax
    Q(:,(j-1)*nB+1:j*nB) = B;
    Z = A(B);
    alpha{j} = B'*Z;
    Z = Z-Q(:,1:j*nB)*(Q(:,1:j*nB)'*Z);
    Z = Z-Q(:,1:j*nB)*(Q(:,1:j*nB)'*Z);
    [B,beta{j}] = qr(Z,0);
end

T = blkdiag(alpha{:});
for ii = 1:j-1
    T((ii-1)*nB+1:ii*nB,ii*nB+1:(ii+1)*nB) = beta{ii}';

    T(ii*nB+1:(ii+1)*nB,(ii-1)*nB+1:ii*nB) = beta{ii};
end
T = (T+T')/2;
Q = Q(:,1:iterMax*nB);

end
