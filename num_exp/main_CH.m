function result = main_CH(name,check,nBmax)
rng(1)
load(name);
A = Problem.A;
n = size(A,2);
eps0 = 1e-3;

D = eps0*rand(n,1);
Af = full(A);
tic
V = null(Af);
result.nulltime = toc;
toc
clear Af V;


result.nBmax = nBmax;
result.eps0 = eps0;
result.check = check;
result.err = zeros(1,nBmax);
result.matvec = zeros(1,nBmax);
result.dim = zeros(1,nBmax);
result.nz = zeros(1,nBmax);
result.time = zeros(1,nBmax);

for nB = nBmax-4:nBmax
b = randn(n,2^(nB-1));
b = orth(b);
tic
[T,Q,result.matvec(nB),result.nz(nB)] = TRlanczos(@(x) A'*(A*x)+D.*x,b,eps0,check);
result.time(nB) = toc;
toc
[V,~]  = eig(T);
result.err(nB) = norm(A*(Q*V(:,1:result.nz(nB))));


end





end