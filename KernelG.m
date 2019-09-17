function[result]=KernelG(x,y,tau)
d=(x-y)*(x-y)';
result=exp(-d/tau^2);
end

