function result = main_GL_LB(nd)


rng(1)
A = prepoA(nd);
n = size(A,1);
eps0 = 1e-3;

D = eps0*rand(n,1);

L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',0.1));

nBmax = 5;

result.eps0 = eps0;
result.matvec = zeros(1,nBmax+1);
result.err = zeros(1,nBmax+1);
result.nz = zeros(1,nBmax+1);
result.time = zeros(1,nBmax+1);
nK = 512;
for ii = 1:nBmax
    nB = 2^(ii+2);
    b = randn(n,nB);
    b = orth(b);
    tic
    [result.matvec(ii),result.nz(ii),result.err(ii)] = TRlanczos2(@(x)A*((L')\x), @(x) L\(A*((L')\x))+D.*x,b,eps0,nK,1.2e-3);
    result.time(ii) = toc;
    toc
end
ii = ii+1;
tic
[result.matvec(ii),result.nz(ii),result.err(ii)] = TRlanczos2(@(x)A*((L')\x), @(x) L\(A*((L')\x)),b,eps0,nK,1.2e-3);
result.time(ii) = toc;
toc

end