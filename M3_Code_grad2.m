clear all
tic
%I1 = double(imread('barbara_contaminated.png'));
I1 = double(imread('cameraman_contaminated.png'));
 [s1 s2 ]=size(I1);

load Omega  %%%%% I1(Omega), I2(Omega) are not contaminated 

% % [Gx,Gy] = imgradientxy(I1,'intermediate');
% Gx1=I1(:,2:256)-I1(:,1:255);
% Gy1=I1(2:256,:)-I1(1:255,:);

IndexD1=zeros(256*255,1);
IndexD2=IndexD1;
IndexD3=linspace(256+1,256*256,256*255)';
IndexD4=linspace(1,256*255,256*255)';
for ii=1:256
    IndexD1(1+(ii-1)*255:255+(ii-1)*255)=linspace(2,256,255)+256*(ii-1);
    IndexD2(1+(ii-1)*255:255+(ii-1)*255)=linspace(1,255,255)+256*(ii-1);
end
    
Index1=linspace(1,256*256,256*256)';

IndexD1=[IndexD3;IndexD1];
IndexD2=[IndexD4;IndexD2];
 
% IndexGx=linspace(1,256*255,256*255);
% IndexGy=linspace(1+256*255,256*255+256*255,256*255);
% Q(IndexGx);

b=A_multi(I1(Index1),Omega);
mu=0.5;
Rho=0.01*mu;
x=I1(Index1);
Q=A_multi(x,IndexD1)-A_multi(x,IndexD2);
V_1=0*b+0.1;
V_2=Q*0;

testerr=zeros(5,3);

for ii=1:5

C=AT_multi(V_1+mu*b,Omega,s1*s2)+AT_multi(V_2+mu*Q,IndexD1,s1*s2)-AT_multi(V_2+mu*Q,IndexD2,s1*s2);
r=mu*AT_multi(x(Omega),Omega,s1*s2)+mu*AT_multi(x(IndexD1),IndexD1,s1*s2)+mu*AT_multi(x(IndexD2),IndexD2,s1*s2)-mu*AT_multi(x(IndexD2),IndexD1,s1*s2)-mu*AT_multi(x(IndexD1),IndexD2,s1*s2)-C;
for jj=1:300
if max(abs(r))>1e-4
alpha=(r'*r)/(r'*(mu*AT_multi(r(Omega),Omega,s1*s2)+mu*AT_multi(r(IndexD1),IndexD1,s1*s2)+mu*AT_multi(r(IndexD2),IndexD2,s1*s2)-mu*AT_multi(r(IndexD2),IndexD1,s1*s2)-mu*AT_multi(r(IndexD1),IndexD2,s1*s2)));
x_new=x-alpha*r;
x=x_new;
r=mu*AT_multi(x(Omega),Omega,s1*s2)+mu*AT_multi(x(IndexD1),IndexD1,s1*s2)+mu*AT_multi(x(IndexD2),IndexD2,s1*s2)-mu*AT_multi(x(IndexD2),IndexD1,s1*s2)-mu*AT_multi(x(IndexD1),IndexD2,s1*s2)-C;
testv=jj;
end
end
testerr(ii,3)=testv;

for kk=1:1 % Q L-2 norm is minimized. 2 iterations will be enough for min QQ.
r=Q+(2/(2+mu))*(V_2-mu*(A_multi(x,IndexD1)-A_multi(x,IndexD2)));
% if max(abs(r))>1e-6
% Q_new=shrink1(x(IndexD1)-x(IndexD2)-V_2/mu,1/mu);
Q_new=Q-r;
Q=Q_new;

% end
max(abs(r));
end

V_1_new=V_1+Rho*(b-x_new(Omega));
V_2_new=V_2+Rho*(Q_new-(x(IndexD1)-x(IndexD2)));

x=x_new;
Q=Q_new;
V_1=V_1_new;
V_2=V_2_new;

testv=A_multi(x,IndexD1)-A_multi(x,IndexD2)-Q;
testerr(ii,1)=sum(abs(testv));
testv=x(Omega)-b;
testerr(ii,2)=sum(abs(testv));
testerr(ii,:);
end
toc
figure;
imshow(reshape(x,s1,s2),[]);