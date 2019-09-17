function[criteria]=findDIM(X,Y,n,p,lambda,numberclose,alphatildeKnorm)
%%It is used to find the useful dimensions%%
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
%%The second derivative%%
M=Secondderivative(n,p,X,Y,w,K);
tuning=zeros(p,1);
for i=1:p
tuning(i)=lambda/alphatildeKnorm(i);
end
max=50;
%%Find the unknown parameter alpha based on given tuning parameter%%
alpha=RKHS(n,p,X,Y,w,K,M,tuning,max);
%%Find the norm for each part of alpha%%
criteria=zeros(p,1);
for i=1:p
criteria(i)=((alpha(((i-1)*n+1):(i*n)))'*(alpha(((i-1)*n+1):(i*n))))^(1/2);
end
end

