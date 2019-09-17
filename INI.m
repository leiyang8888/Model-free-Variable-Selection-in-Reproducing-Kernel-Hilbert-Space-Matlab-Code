function [K,M,w,alphatildeKnorm] = INI(X,Y,n,p,numberclose)
%%Find the alphatilde to give the adaptive lasso weight%%
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
tuning=zeros(p,1)+0.000001;
%%50 iteration to get estimation for alphatilde%%
betatilde=RKHS(n,p,X,Y,w,K,M,tuning,50);
alphatildeKnorm=zeros(p,1);
%%Find the corresponding K norm%%
for i=1:p
alphatildeKnorm(i)=betatilde(((i-1)*n+1):(i*n))'*K*betatilde(((i-1)*n+1):(i*n));
end
end

