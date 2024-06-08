function A = prepoA(nd)
A = mmread('graph.mtx');
idx = sum(A,1)>=nd;
A = A(idx,idx);
D = sum(A,1);
idx = find(D>0);
A = A(idx,idx);
D = sum(A,1);
A = diag(D)-A;
end