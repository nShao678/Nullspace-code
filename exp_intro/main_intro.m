clear all
rng(1)
m = 64;
n = 10000+m;
A = [zeros(1,m),1+(1:n-m)/n]';
X = randn(n,m);
blk = [1,2,4,8,16,32,64];
len = length(blk);
nz = zeros(1,len);
iterMax = 1024;
for ii = 1:len
    X0 = orth(X(:,1:blk(ii)));
    [T,Q,nz(ii)] = Blanczos(@(x) A.*x,X0,iterMax);
end
latex(sym([blk;nz]))

