clear all
tic
I1 = double(imread('barbara_contaminated.png')); % I1 is a 256 by 256 matrix
% I1 = double(imread('cameraman_contaminated.png'));
load Omega  %%%%% I1(Omega), I2(Omega) are not contaminated. Omega is the index of uncontaminated pixels
level=2;% parameter of wavelet transformation function

%%%%%%   Wavelet transformation  W*u,, W^T*W = I
coef = swt2(I1,1,level);% this is Wx
[s1 s2 s3]=size(coef);
%%% Soft Thresholding needs to apply on coef(:,:,2:end)
%%%%%  Inverse wavelet transformation W^T*u
% newI = iswt2(coef,1,1);

Index1=linspace(1,s1*s2,s1*s2)';
Index2=linspace(1,s1*s2*s3,s1*s2*s3)';

I1=I1(Index1);%chenge matrix to vector
% Q=coef(Index2);
b=I1(Omega);%calculate one constraint(Ax=b)
mu=0.05;% parameter in ALM
Rho=0.01*mu;% step size in ALM with Rho < mu
x=I1*0;% initial condiction of x
V_1=0*b+0.1;% initial condiction of V_1
V_2=0*coef(Index2);% initial condiction of V_2
temp=swt2(reshape(x,s1,s2),1,level);% calculate Ax=Q(Q is in matrix)
Q=temp(Index2);% change Q to a vector form
testerr=zeros(500,3);% calculation errors and iteration numbers

for ii=1:500 %number of iteration

C1=AT_multi(V_1+mu*b,Omega,s1*s2);
C2=iswt2(reshape(V_2+mu*Q,s1,s2,s3),1,level);
C2=C2(Index1);
C=C1+C2; %C is calculated

% Newton mathod to  min x(k+1) 
r=mu*AT_multi(x(Omega),Omega,s1*s2)+mu*x-C; % r=Hx-c
for jj=1:50 % max number of iteration
if max(abs(r))>1e-8 % error 
alpha=(r'*r)/(r'*(mu*AT_multi(r(Omega),Omega,s1*s2)+mu*r)); %alpha=rr/rHr
x_new=x-alpha*r;
x=x_new;
r=mu*AT_multi(x(Omega),Omega,s1*s2)+mu*x-C; % update r
testv=jj;
end
end
testerr(ii,3)=testv;

temp=swt2(reshape(x_new,s1,s2),1,level);
temp=temp(Index2);
Q_new=shrink1(temp-V_2/mu,1/mu);% min Q

V_1_new=V_1+Rho*(b-x_new(Omega)); % calculate new V_1
V_2_new=V_2+Rho*(Q_new-temp); % calculate new V_2

x=x_new;
Q=Q_new;
V_1=V_1_new;
V_2=V_2_new;% update x Q V_1 V_2

testv=swt2(reshape(x,s1,s2),1,level);
testv=testv(Index2);
testv=Q-testv; %calculate error of Q-wx
testerr(ii,1)=sum(abs(testv));
testv=x(Omega)-b;
testerr(ii,2)=sum(abs(testv)); %calculate error of Ax-b
testerr(ii,:);
end
toc 
figure;
imshow(reshape(x,s1,s2),[]);


