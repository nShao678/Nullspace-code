function result = main_GL_SV(eps0,nd)


rng(1)
A = prepoA(nd);
n = size(A,1);


D = eps0*rand(n,1);

b = randn(n,1);
b = orth(b);
L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',0.1));
iterMax = 1000;

[T,Q] = lanczos(@(x) L\(A*((L')\x))+D.*x,b,iterMax);




result.eps0 = eps0;
[V,eval] = eig(T,'vector');
nz = sum(eval<4*eps0);
result.nz = nz;
result.lambdaN = eval(nz);
Ql = L'\Q;
V = Ql*V(:,1:nz);
M = V'*(A*V);
M = (M+M')/2;
result.normB = norm(M);
end