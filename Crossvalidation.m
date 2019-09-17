function[STAB,CRITERIA1,CRITERIA2]=Crossvalidation(X,Y,n,p,n1,n2,numberclose,totaltime,alphatildeKnorm)
%%It is used to choose the tuning parameter based on training data%%
t=20;
STAB=zeros(totaltime,t);
CRITERIA1=zeros(p,t,totaltime);
CRITERIA2=zeros(p,t,totaltime);
%%The step length and start point for tuning parameter%%
start=-2;
step=4;
%%Do the cross validation for totaltime times%%
for time=1:totaltime
%%Randomly divide the data into two parts%%
X1=zeros(n1,p);
Y1=zeros(n1,1);
X2=zeros(n2,p);
Y2=zeros(n2,1);
xulie=randperm(n);
for j=1:n1
X1(j,:)=X(xulie(j),:);
Y1(j)=Y(xulie(j));
end
for j=1:n2
X2(j,:)=X(xulie(j+n1),:);
Y2(j)=Y(xulie(j+n1));
end
K1=zeros(n1,n1);
%%Find the Gaussian Kernel matrix%%
Distance=zeros(n1,n1);
for i=1:n1
for j=1:n1
Distance(i,j)=sqrt((X1(i,:)-X1(j,:))*(X1(i,:)-X1(j,:))');
end
end
tau1=median(median(Distance));
for i=1:n1
for j=1:n1
K1(i,j)=KernelG(X1(i,:),X1(j,:),tau1);
end
end
K2=zeros(n2,n2);
Distance=zeros(n2,n2);
for i=1:n2
for j=1:n2
Distance(i,j)=sqrt((X2(i,:)-X2(j,:))*(X2(i,:)-X2(j,:))');
end
end
tau2=median(median(Distance));
for i=1:n2
for j=1:n2
K2(i,j)=KernelG(X2(i,:),X2(j,:),tau2);
end
end
%%Compute the weight and decide how many close points will be selected for these two data sets%%
w1=Findweight(n1,p,X1,Y1,tau1,numberclose);
M1=Secondderivative(n1,p,X1,Y1,w1,K1);
w2=Findweight(n2,p,X2,Y2,tau2,numberclose);
M2=Secondderivative(n2,p,X2,Y2,w2,K2);
max=10;
tuning1=zeros(p,1);
tuning2=zeros(p,1);
for t0=1:t
VS1=zeros(p,1);
VS2=zeros(p,1);
lambda=10^(start+(t0-1)*step/19);
%%The adaptive tuning parameter%%
for i=1:p
tuning1(i)=lambda/alphatildeKnorm(i);
tuning2(i)=lambda/alphatildeKnorm(i);
end
%%Get the estimation for alpha for these two data sets%%
alpha1=RKHS(n1,p,X1,Y1,w1,K1,M1,tuning1,max);
alpha2=RKHS(n2,p,X2,Y2,w2,K2,M2,tuning2,max);
%%The dimension choose for each data set%%
criteria1=zeros(p,1);
criteria2=zeros(p,1);
for i=1:p
criteria1(i)=((alpha1(((i-1)*n1+1):(i*n1)))'*(alpha1(((i-1)*n1+1):(i*n1))))^(1/2);
criteria2(i)=((alpha2(((i-1)*n2+1):(i*n2)))'*(alpha2(((i-1)*n2+1):(i*n2))))^(1/2);
if criteria1(i)>=0.001
VS1(i)=1;
else
end
if criteria2(i)>=0.001
VS2(i)=1;
else
end
end
n11=0;
n12=0;
n21=0;
n22=0;
for i=1:p
if VS1(i)==1 && VS2(i)==1
n11=n11+1;
else
end
if VS1(i)==1 && VS2(i)==0
n12=n12+1;
else
end
if VS1(i)==0 && VS2(i)==1
n21=n21+1;
else
end
if VS1(i)==0 && VS2(i)==0
n22=n22+1;
else
end
end
%%Find the variable selection stability based on these two data sets%%
Pra=(n11+n22)/p;
Pre=(n11+n12)*(n11+n21)/p^2+(n12+n22)*(n21+n22)/p^2;
STAB(time,t0)=(Pra-Pre)/(1-Pre);
if sum(VS1)==p
STAB(time,t0)=-1;   
else    
end
if sum(VS2)==p
STAB(time,t0)=-1;   
else    
end
if sum(VS1)==0
STAB(time,t0)=-1;   
else    
end
if sum(VS2)==0
STAB(time,t0)=-1;   
else    
end
CRITERIA1(:,t0,time)=criteria1;
CRITERIA2(:,t0,time)=criteria2;
end
end
end

