function[alpha,iteration]=RKHS(n,p,X,Y,w,K,M,tuning,max)
%%It is used to renew the estimation for unknown parameter%%
%%Find the largest eigenvalue for each block of the second derivative%%
gamma=largesteigenvalue(M,n,p);
alpha=zeros(n*p,1);
eps1=10;
iteration=0;
%%Max is the maximum number of iteration%%
while(eps1>10^(-6) && iteration<=max)
b1=alpha;
%%Renew alpha based on the formula%%
for l=1:p
alphal=alpha(((l-1)*n+1):(l*n));
delta=alphal*gamma(l)-Firstderivative(n,p,X,Y,w,K,alpha,M,l);
partone=delta/gamma(l);
parttwo=1-tuning(l)/(delta'*delta)^(1/2);
if(parttwo<0) 
parttwo=0;
else
end
alphanew=parttwo*partone;
alpha(((l-1)*n+1):(l*n))=alphanew;
end
b2=alpha;
%%The criteria for stopping renew process%%
eps1=((b1-b2)'*(b1-b2))^(1/2);
iteration=iteration+1;
end
end

