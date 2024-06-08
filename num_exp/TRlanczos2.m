function [count,nz,err] = TRlanczos2(A0,A,B,eps0,nK,tol)

% initial with lanczos
[n,nB] = size(B);
iternK = nK/nB;
B = orth(B);

epstol = sqrt(eps);
epsmin = eps*nB*sqrt(n);

numRO = 0;
count = 0;
keep = 0;
Theta = [];
Y = [];
W = [];
nz = 0;
for ii = 1:n
    alpha = cell((nK-keep)/nB,1);
    beta = cell((nK-keep)/nB,1);
    iterkeep = keep/nB;
    Q = [Y,B];
    Z = A(B);
    count = count+nB;
    tmp = B'*Z;
    tmp = (tmp+tmp')/2;
    alpha{1} = [diag(Theta),W';W,tmp];
    Z = Z-Q*(Q'*Z);
    % Z = Z-Q*(Q'*Z);
    [B,beta{1}] = qr(Z,0);
    Q = [Q,B,zeros(n,nK-keep-2*nB)];
    normalpha = zeros(1,iternK);
    normbeta = zeros(1,iternK);
    for j = 1:iterkeep
        idx = (j-1)*nB+1:j*nB;
        normalpha(j) = max(abs(Theta(idx)));
        normbeta(j) = norm(W(:,idx));
    end
    normalpha(iterkeep+1) = norm(tmp);
    normbeta(iterkeep+1) = norm(beta{1});
    omega = zeros(iternK+1,iternK);
    omega(2:iterkeep+2,1:iterkeep+1) = tril(epsmin*ones(iterkeep+1));
    %lanczos
    for j = 2:(nK-keep)/nB
        jj = iterkeep+j;
        Z = A(B);
        count = count+nB;
        Z = Z-Q(:,(jj-2)*nB+1:(jj-1)*nB)*beta{j-1}';
        alpha{j} = B'*Z;
        alpha{j} = (alpha{j}+alpha{j}')/2;
        Z = Z-B*alpha{j};
%         tmp = B'*Z;
%         tmp = (tmp+tmp')/2;
%         alpha{j} = alpha{j}+tmp;
%         Z = Z-B*tmp;
        normalpha(jj) = norm(alpha{j});
        [B,beta{j}] = qr(Z,0);
        normbeta(jj) = norm(beta{j});
        Q(:,jj*nB+1:(jj+1)*nB) = B;


        ninvj = normi(beta{j});
        omega(jj+1,jj) = epsmin;
        omega(jj+1,iterkeep+1) = (normalpha(jj)+normalpha(iterkeep+1))*omega(jj,iterkeep+1)+normbeta(iterkeep+1)*omega(jj,iterkeep+2)+normbeta(jj-1)*omega(jj-1,iterkeep+1)+2*epsmin;
        for l = 1:iterkeep
            omega(jj+1,l) = ninvj*(normbeta(l)*omega(l,iterkeep+1)+normbeta(jj-1)*omega(jj-1,l)+(normalpha(jj)+normalpha(l))*omega(jj,l)+2*epsmin);
            omega(jj+1,iterkeep+1) = omega(jj+1,iterkeep+1)+normbeta(l)*omega(jj,l);
        end
        omega(jj+1,iterkeep+1) = ninvj*omega(jj+1,iterkeep+1);
        for l = iterkeep+2:jj-1
            omega(jj+1,l) = ninvj*(normbeta(l)*omega(jj,l+1)+normbeta(l-1)*omega(jj,l-1)+normbeta(jj-1)*omega(jj-1,l)+(normalpha(jj)+normalpha(l))*omega(jj,l)+2*epsmin);
        end
        if max(omega(jj+1,:))>epstol
            Q(:,(jj-1)*nB+1:jj*nB) = Q(:,(jj-1)*nB+1:jj*nB)-Q(:,1:(jj-1)*nB)*(Q(:,1:(jj-1)*nB)'*Q(:,(jj-1)*nB+1:jj*nB));
%             Q(:,(jj-1)*nB+1:jj*nB) = Q(:,(jj-1)*nB+1:jj*nB)-Q(:,1:(jj-1)*nB)*(Q(:,1:(jj-1)*nB)'*Q(:,(jj-1)*nB+1:jj*nB));
            [Q0,R0] = qr(Q(:,(jj-1)*nB+1:jj*nB),0);
            Q(:,(jj-1)*nB+1:jj*nB) = Q0.*sign(diag(R0))';
%             beta{j-1} = R0*beta{j-1};
%             normbeta(jj-1) = norm(beta{j-1});
            Z = Z-Q(:,1:jj*nB)*(Q(:,1:jj*nB)'*Z);
%             Z = Z-Q(:,1:jj*nB)*(Q(:,1:jj*nB)'*Z);
            [B,beta{j}] = qr(Z,0);
            normbeta(jj) = norm(beta{j});
            Q(:,jj*nB+1:(jj+1)*nB) = B;
            numRO = numRO+1;
            omega(jj,1:jj-1) = epsmin;
            omega(jj+1,1:jj) = epsmin;
        end
    end
    T = blkdiag(alpha{:});
    for l = 1:j-1
        T((l-1)*nB+keep+1:l*nB+keep,l*nB+keep+1:(l+1)*nB+keep) = beta{l}';
        T(l*nB+keep+1:(l+1)*nB+keep,(l-1)*nB+keep+1:l*nB+keep) = beta{l};
    end
    [S,Theta] = eig(T,'vector');
    nz = sum(Theta<3*eps0);
    err = norm(A0(Q(:,1:nK)*S(:,1:nz)));
    if nz>=1 && err<tol
        break
    end
    keep = floor((nz+nK)/2/nB)*nB;
    keep = min(keep,nK-3*nB);
    Y = Q(:,1:nK)*S(:,1:keep);

    W = zeros(nB,nK);
    W(:,nK-nB+1:nK) = beta{j};
    W = W*S(:,1:keep);

    B = B-Y*(Y'*B);
    % B = B-Y*(Y'*B);
    Theta = Theta(1:keep);
    
end
end
function y = normi(A)
y = 1/min(svd(A));
end