clear all
m = 21;
n = 420;
A = [zeros(1,m),(1:n-m)]';
X = randn(n,m);
blk = [1,2,3,4,5,6,7,10,12,14,15,20,21];
len = length(blk);
nz = zeros(1,len);
iterMax = n;
for ii = 1:len
    X0 = orth(X(:,1:blk(ii)));
    [~,~,nz(ii)] = Blanczos(@(x) A.*x,X0,iterMax);
end
latex(sym([blk;nz]))

