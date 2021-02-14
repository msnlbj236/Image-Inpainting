%function ReadData

I1 = double(imread('barbara_contaminated.png'));
I2 = double(imread('cameraman_contaminated.png'));
figure;
imshow(I1,[]);
figure;
imshow(I2,[]);

load Omega  %%%%% I1(Omega), I2(Omega) are not contaminated 



%%%%%%   Wavelet transformation  W*u,, W^T*W = I
coef = swt2(I2,1,1);
%%% Soft Thresholding needs to apply on coef(:,:,2:end)
%%%%%  Inverse wavelet transformation W^T*u
newI = iswt2(coef,1,1);

figure;
imshow(newI,[]);
