function result = main_GL_SVH(eps0,nd)


rng(1)
A = prepoA(nd);
n = size(A,1);


D = eps0*rand(n,1);

b = randn(n,1);
b = orth(b);
L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',0.1));

iterMax = 1000;

[T,~] = lanczos(@(x) L\(A*((L')\x))+D.*x,b,iterMax);




result.eps0 = eps0;
eval = eig(T);
nz = sum(eval<3*eps0);
result.nz = nz;
jj = 1;
result.zero = zeros(1,nz);
for ii = 1:size(T,1)
    eval = eig(T(1:ii,1:ii));
    nzt = sum(eval<3*eps0);
    if nzt==jj
        result.zero(jj) = ii;
        jj = jj+1;
    end
end
end