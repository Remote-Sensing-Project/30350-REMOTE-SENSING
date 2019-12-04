clc
close all
%%
% Script for plotting The Backscatter and correlation coefficient

matrixB = get_data;
DiaB = zeros([3,1048576]);
HH = zeros([1,1048576]);
VV = zeros([1,1048576]);
HV = zeros([1,1048576]);
corr = zeros([1,1048576]);
% A is HH, B is VV ,C is HV and corr is HHVV matix solved. 
for n = 1:1048576
    DiaB(:,n)  = diag(matrixB(:,:,n));
    HH(:,n) = DiaB(1,n);
    VV(:,n) = DiaB(2,n);
    HV(:,n) = DiaB(3,n);
    HHVV(:,n) = matrixB(1,3,n);
    corr(:,n) = HHVV(:,n)/(HH(:,n)*VV(:,n))^(1/2); 
end
%%
% ---Reshaping the Data---- %
ImgA = reshape(HH,1024,1024);   % reshaping HH  
ImgB = reshape(VV,1024,1024);   % reshaping VV
ImgC = reshape(HV,1024,1024);   % reshaping HV
ImgD = reshape(abs(corr),1024,1024); % reshaping the basolute value of CC
%%
%-------------------------------------------------------------------------%
subplot(2,2,1)      %-- Sub plott 1 --%
imagesc(log10(ImgA))% Plotting HH
title('HH');        % Title

subplot(2,2,2)      %-- Sub plott 2 --%
imagesc(log10(ImgB))% Plottting VV
title('VV');        % Title

subplot(2,2,3)      %-- Sub plott 3 --%
imagesc(log10(ImgC))% Plotting HV
title('HV');        % Title

subplot(2,2,4)      %-- Sub plott 4 --%
imagesc(abs(ImgD))  % Plotting correlation coefficient
title('correlation coefficient'); % Title