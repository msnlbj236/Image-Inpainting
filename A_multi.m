function [outputArg1] = A_multi(x,Omega)
%A_MULTI 

N=zeros(length(Omega),1);

for ii=1:length(Omega)

N(ii)=x(Omega(ii));
end
outputArg1=N;
end

