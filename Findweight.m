function[w]=Findweight(n,p,X,Y,tau,numberclose)
%%It is used to find the weight for each observation%%
w=zeros(n,n);
weight=w;
%%Find the distance%%
for i=1:n
d=zeros(n,1);
for j=1:n
d(j)=(X(i,:)-X(j,:))*(X(i,:)-X(j,:))';
if i==j
d(j)=Inf;
else
end
end
%%Choose how many close data points%%
d1=sort(d);
cutoff=d1(numberclose);
for j=1:n
if d(j)<=cutoff
weight(i,j)=exp(-d(j)/tau^2);
else
end
end
end
w=weight;
end

