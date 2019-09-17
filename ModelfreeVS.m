function [ criteria,STAB ] = ModelfreeVS( X,Y,numberclose,totaltime)
%%The function here is used to find the informative random variables.
%%Input: (1) X:variables data set (2)Y:response (3)numberclose:how
%%many around data points will be introduced into weighted loss function
%%(4)totaltime:The times for doing crossvalidation
[n p]=size(X);
n1=floor(n/2);
n2=floor(n/2);
%%Find the adaptive weight%%
[K,M,w,alphatildeKnorm]=INI(X,Y,n,p,numberclose);
%%Variable selection stability for a sequence of tuning parameters%%
[STAB,CRITERIA1,CRITERIA2]=Crossvalidation(X,Y,n,p,n1,n2,numberclose,totaltime,alphatildeKnorm);

%%Find the most reasonable tuning parameter%%
if max(mean(STAB))>0
biaozhun=max(mean(STAB))*0.9;
tt=find(mean(STAB)>=biaozhun);
t=tt(1);
else
end
if max(mean(STAB))<=0
biaozhun=max(mean(STAB));
tt=find(mean(STAB)==biaozhun);
t=tt(1);
else
end

%%Decide which variable is informative%%
lambda=10^(-2+(t-1)*4/19);
X3=X(1:floor(n/2),:);
Y3=Y(1:floor(n/2),:);
[criteria]=findDIM(X3,Y3,n/2,p,lambda,numberclose,alphatildeKnorm);

end

