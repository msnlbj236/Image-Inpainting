function [y] = shrink1(x,t)
%SHRINK1 
M=0*x;
for ii=1:length(x)
    if x(ii)>t
        M(ii)=x(ii)-t;
    end
    if x(ii)<-1*t
        M(ii)=x(ii)+t;
    end
     
end
y=M;
