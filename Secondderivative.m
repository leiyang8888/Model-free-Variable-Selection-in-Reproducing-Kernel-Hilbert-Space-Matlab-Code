function[M]=Secondderivative(n,p,X,Y,w,K)
%%Find the second derivative for the l-th block component of alpha%%
M=zeros(n*p,n*p);
for i=1:n
for j=1:n
if w(i,j)>0;
A=(X(i,:)-X(j,:))'*(X(i,:)-X(j,:));
B=K(:,i)*K(:,i)';
M=M+2*w(i,j)*kron(A,B)/(n*(n-1));
else
end
end
end
end

