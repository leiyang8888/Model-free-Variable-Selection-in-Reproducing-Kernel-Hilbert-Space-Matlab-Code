function[gamma]=largesteigenvalue(M,n,p)
%%Find the largest eigenvalue of the n*n block diagonal matrix of M%%
gamma=zeros(p,1);
A=zeros(n,n);
for i=1:p
A=M(((i-1)*n+1):(i*n),((i-1)*n+1):(i*n));
gamma(i)=max(eig(A));
end
end

