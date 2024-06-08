function result = main_CH_LB
rng(1)
load('SSmatrix\Franz9.mat');
A = Problem.A;
n = size(A,2);
eps0 = 1e-3;

D = eps0*rand(n,1);

iiMax = 5;
result.eps0 = eps0;
result.err = zeros(1,iiMax+1);
result.matvec = zeros(1,iiMax+1);
result.nz = zeros(1,iiMax+1);
result.time = zeros(1,iiMax+1);

nK = 2560;
for ii = 1:iiMax
    nB = 2^(ii)*20;
    b = randn(n,nB);
    b = orth(b);
    
    tic
    [result.matvec(ii),result.nz(ii),result.err(ii)] = TRlanczos2(@(x) A*x,@(x) A'*(A*x)+D.*x,b,eps0,nK,1e-4);
    result.time(ii) = toc;
    toc
end
ii = ii+1;
tic
[result.matvec(ii),result.nz(ii),result.err(ii)] = TRlanczos2(@(x) A*x,@(x) A'*(A*x),b,eps0,nK,1e-4);
result.time(ii) = toc;
toc


end