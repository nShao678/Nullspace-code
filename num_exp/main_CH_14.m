function result = main_CH_14
rng(1)
load('SSmatrix\GL7d14.mat');
A = Problem.A;
n = size(A,2);
eps0 = 1e-3;



D = eps0*rand(n,1);

iiMax = 6;
result.eps0 = eps0;
result.err = zeros(1,iiMax);
result.matvec = zeros(1,iiMax);
result.nz = zeros(1,iiMax);
result.time = zeros(1,iiMax);
nK = 10240;
for ii = 1:iiMax
    nB = 2^(ii+3)
    b = randn(n,nB);
    b = orth(b);
    tic
    [T,Q,result.matvec(ii),result.nz(ii)] = TRlanczos(@(x) A'*(A*x)+D.*x,b,eps0,nK);
    result.time(ii) = toc;
    toc
    [V,~]  = eig(T);
    result.err(ii) = norm(A*(Q*V(:,1:result.nz(ii))));
end


end