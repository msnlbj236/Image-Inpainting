function [outputArg1] = AT_multi(x,Omega,num)
%A_MULTI 

N=zeros(num,1);

for ii=1:length(Omega)

N(Omega(ii))=x(ii)+N(Omega(ii));
end

outputArg1=N;
end

