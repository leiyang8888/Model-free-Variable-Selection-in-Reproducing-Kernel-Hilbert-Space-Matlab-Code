function[Attribute,Response,NORM]=INI_ridge(X,Y,numberclose)
%%It is used to find the adaptive lasso weight via ridge regression%%
[n p]=size(X);
K=zeros(n,n);
%%Find the Gaussian Kernel matrix%%
Distance=zeros(n,n);
for i=1:n
for j=1:n
Distance(i,j)=sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
end
end
tau=median(median(Distance));
for i=1:n
for j=1:n
K(i,j)=KernelG(X(i,:),X(j,:),tau);
end
end
%%Find the weight for each training data point%%
w=Findweight(n,p,X,Y,tau,numberclose);
Attribute=zeros(n*p,n*p);
Response=zeros(n*p,1);

%%Find the response and design matrix for ridge regression%%
for i=1:n
for j=1:n
if w(i,j)>0
Attribute=Attribute+w(i,j)*(kron(X(i,:)-X(j,:),K(i,:)))'*kron(X(i,:)-X(j,:),K(i,:))/n/(n-1);
Response=Response+w(i,j)*(Y(i)-Y(j))*(kron(X(i,:)-X(j,:),K(i,:)))'/n/(n-1);
else
end
end
end

lambda=10^(-6);
result=pinv(Attribute+lambda*kron(eye(p),K))*Response;
NORM=zeros(p,1);
for i=1:p
NORM(i)=result((i-1)*n+1:i*n)'*K*result((i-1)*n+1:i*n);
end

end

