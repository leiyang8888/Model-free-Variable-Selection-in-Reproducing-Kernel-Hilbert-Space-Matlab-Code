function[result]=Firstderivative(n,p,X,Y,w,K,alpha,M,l)
%%Find the first derivative for the l-th block component of alpha%%
result=zeros(n,1);
U=zeros(n,1);
for i=1:n
for j=1:n
if w(i,j)>0
U=U+2*w(i,j)*(Y(i)-Y(j))*(X(i,l)-X(j,l))*K(i,:)'/(n*(n-1));
else
end
end
end
result=result+M(((l-1)*n+1):(l*n),:)*alpha-U;
end

