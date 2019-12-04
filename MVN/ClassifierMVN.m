%
% Wishart Classification script:
% Multi-frequency classifyer for polarometric data
% 
clc
%clear all
close all

% Data set 1
matrixB1 = get_data_1();

% Reference data for the creation of the masks
DiaB = zeros([3,1048576]);
HH = zeros([1,1048576]);
VV = zeros([1,1048576]);
HV = zeros([1,1048576]);
corr = zeros([1,1048576]);
for n = 1:1048576
    DiaB(:,n)  = diag(matrixB1(:,:,n));
    HH(:,n) = DiaB(1,n);
    HV(:,n) = DiaB(2,n);
    VV(:,n) = DiaB(3,n);
    vP(:,n) = HH(:,n) + VV(:,n) + 2*HV(:,n);
end

% Data set 2
matrixB2 = get_data_2();
% Data set 3
matrixB3 = get_data_3();
% Data set 4
matrixB4 = get_data_4();
% Data set 5
matrixB5 = get_data_5();
% Data set 6
matrixB6 = get_data_6();
% Data set 7
matrixB7 = get_data_7();
% Data set 8
matrixB8 = get_data_8();
%% Processing the data (10 classes and backscatter coeff's) 
%Getting the masks
Nclasses = 10;
% masks = zeros([1024,1024,Nclasses]);

for k = 1:Nclasses

    %Getting the mask
%     figure 
%     ImgvP = reshape(vP,1024,1024);
%     masks(:,:,k) = roipoly(ImgvP);
    q = find(masks(:,:,k) == 1);

    prmeanU1 = zeros([3,1]);
    rMatrixB1 = matrixB1(:,:,q);
    % sigma1 = zeros([3, 3]);

    prmeanU2 = zeros([3,1]);
    rMatrixB2 = matrixB2(:,:,q);
    % sigma2 = zeros([3, 3]);

    prmeanU3 = zeros([3,1]);
    rMatrixB3 = matrixB3(:,:,q);
    % sigma3 = zeros([3, 3]);

    prmeanU4 = zeros([3,1]);
    rMatrixB4 = matrixB4(:,:,q);
    %sigma4 = zeros([3, 3]);

    prmeanU5 = zeros([3,1]);
    rMatrixB5 = matrixB5(:,:,q);
    % sigma5 = zeros([3, 3]);
    
    prmeanU6 = zeros([3,1]);
    rMatrixB6 = matrixB6(:,:,q);
    % sigma6 = zeros([3, 3]);
    
    prmeanU7 = zeros([3,1]);
    rMatrixB7 = matrixB7(:,:,q);
    % sigma7 = zeros([3, 3]);c6dd
    
    prmeanU8 = zeros([3,1]);
    rMatrixB8 = matrixB8(:,:,q);
    % sigma8 = zeros([3, 3]);


    %Getting u and sigma in the mask
    for n = 1:length(q)
        DiaB1M(:,n)  = diag(rMatrixB1(:,:,n));
        HH1M(:,n) = DiaB1M(1,n);
        HV1M(:,n) = DiaB1M(2,n);
        VV1M(:,n) = DiaB1M(3,n);
        u1M(:,n) = [log10(HH1M(:,n)); log10(HV1M(:,n)); log10(VV1M(:,n))];
        prmeanU1 = prmeanU1 + u1M(:,n);

        DiaB2M(:,n)  = diag(rMatrixB2(:,:,n));
        HH2M(:,n) = DiaB2M(1,n);
        HV2M(:,n) = DiaB2M(2,n);
        VV2M(:,n) = DiaB2M(3,n);
        u2M(:,n) = [log10(HH2M(:,n)); log10(HV2M(:,n)); log10(VV2M(:,n))];
        prmeanU2 = prmeanU2 + u2M(:,n);

        DiaB3M(:,n)  = diag(rMatrixB3(:,:,n));
        HH3M(:,n) = DiaB3M(1,n);
        HV3M(:,n) = DiaB3M(2,n);

        VV3M(:,n) = DiaB3M(3,n);
        u3M(:,n) = [log10(HH3M(:,n)); log10(HV3M(:,n)); log10(VV3M(:,n))];
        prmeanU3 = prmeanU3 + u3M(:,n);

        DiaB4M(:,n)  = diag(rMatrixB4(:,:,n));
        HH4M(:,n) = DiaB4M(1,n);
        HV4M(:,n) = DiaB4M(2,n);
        VV4M(:,n) = DiaB4M(3,n);
        u4M(:,n) = [log10(HH4M(:,n)); log10(HV4M(:,n)); log10(VV4M(:,n))];
        prmeanU4 = prmeanU4 + u4M(:,n);

        DiaB5M(:,n)  = diag(rMatrixB5(:,:,n));
        HH5M(:,n) = DiaB5M(1,n);
        HV5M(:,n) = DiaB5M(2,n);
        VV5M(:,n) = DiaB5M(3,n);
        u5M(:,n) = [log10(HH5M(:,n)); log10(HV5M(:,n)); log10(VV5M(:,n))];
        prmeanU5 = prmeanU5 + u5M(:,n);
    
        DiaB6M(:,n)  = diag(rMatrixB6(:,:,n));
        HH6M(:,n) = DiaB6M(1,n);
        HV6M(:,n) = DiaB6M(2,n);
        VV6M(:,n) = DiaB6M(3,n);
        u6M(:,n) = [log10(HH6M(:,n)); log10(HV6M(:,n)); log10(VV6M(:,n))];
        prmeanU6 = prmeanU6 + u6M(:,n);
    
        DiaB7M(:,n)  = diag(rMatrixB7(:,:,n));
        HH7M(:,n) = DiaB7M(1,n);
        HV7M(:,n) = DiaB7M(2,n);
        VV7M(:,n) = DiaB7M(3,n);
        u7M(:,n) = [log10(HH7M(:,n)); log10(HV7M(:,n)); log10(VV7M(:,n))];
        prmeanU7 = prmeanU7 + u7M(:,n);
    
        DiaB8M(:,n)  = diag(rMatrixB8(:,:,n));
        HH8M(:,n) = DiaB8M(1,n);
        HV8M(:,n) = DiaB8M(2,n);
        VV8M(:,n) = DiaB8M(3,n);
        u8M(:,n) = [log10(HH8M(:,n)); log10(HV8M(:,n)); log10(VV8M(:,n))];
        prmeanU8 = prmeanU8 + u8M(:,n);
    end

    meanU1 = (1/length(q))*(prmeanU1);
    meanU2 = (1/length(q))*(prmeanU2);
    meanU3 = (1/length(q))*(prmeanU3);
    meanU4 = (1/length(q))*(prmeanU4);
    meanU5 = (1/length(q))*(prmeanU5);
    meanU6 = (1/length(q))*(prmeanU6);
    meanU7 = (1/length(q))*(prmeanU7);
    meanU8 = (1/length(q))*(prmeanU8);

    meanU = meanU1;
    meanU(4:6) = meanU2;
    meanU(7:9) = meanU3;
    meanU(10:12) = meanU4;
    meanU(13:15) = meanU5;
    meanU(16:18) = meanU6;
    meanU(19:21) = meanU7;
    meanU(22:24) = meanU8;


    % Calculating the pixel probabilities

    for n = 1:1048576

        %WE SHOULD BE ABLE TO GET ALL DATA SEPARATELY TO MAKE THIS MORE EFFICIENT(DIFFERENT SCRIPT)
        DiaB1(:,n)  = diag(matrixB1(:,:,n));
        HH1(:,n) = DiaB1(1,n);
        HV1(:,n) = DiaB1(2,n);
        VV1(:,n) = DiaB1(3,n);
        u1(:,n) = [log10(HH1(:,n)); log10(HV1(:,n)); log10(VV1(:,n))];

        DiaB2(:,n)  = diag(matrixB2(:,:,n));
        HH2(:,n) = DiaB2(1,n);
        HV2(:,n) = DiaB2(2,n);
        VV2(:,n) = DiaB2(3,n);
        u2(:,n) = [log10(HH2(:,n)); log10(HV2(:,n)); log10(VV2(:,n))];

        DiaB3(:,n)  = diag(matrixB3(:,:,n));
        HH3(:,n) = DiaB3(1,n);
        HV3(:,n) = DiaB3(2,n);
        VV3(:,n) = DiaB3(3,n);
        u3(:,n) = [log10(HH3(:,n)); log10(HV3(:,n)); log10(VV3(:,n))];

        DiaB4(:,n)  = diag(matrixB4(:,:,n));
        HH4(:,n) = DiaB4(1,n);
        HV4(:,n) = DiaB4(2,n);
        VV4(:,n) = DiaB4(3,n);
        u4(:,n) = [log10(HH4(:,n)); log10(HV4(:,n)); log10(VV4(:,n))];

        DiaB5(:,n)  = diag(matrixB5(:,:,n));
        HH5(:,n) = DiaB5(1,n);
        HV5(:,n) = DiaB5(2,n);
        VV5(:,n) = DiaB5(3,n);
        u5(:,n) = [log10(HH5(:,n)); log10(HV5(:,n)); log10(VV5(:,n))];
        
        DiaB6(:,n)  = diag(matrixB6(:,:,n));
        HH6(:,n) = DiaB6(1,n);
        HV6(:,n) = DiaB6(2,n);
        VV6(:,n) = DiaB6(3,n);
        u6(:,n) = [log10(HH6(:,n)); log10(HV6(:,n)); log10(VV6(:,n))];
        
        DiaB7(:,n)  = diag(matrixB7(:,:,n));
        HH7(:,n) = DiaB7(1,n);
        HV7(:,n) = DiaB7(2,n);
        VV7(:,n) = DiaB7(3,n);
        u7(:,n) = [log10(HH7(:,n)); log10(HV7(:,n)); log10(VV7(:,n))];
        
        DiaB8(:,n)  = diag(matrixB8(:,:,n));
        HH8(:,n) = DiaB8(1,n);
        HV8(:,n) = DiaB8(2,n);
        VV8(:,n) = DiaB8(3,n);
        u8(:,n) = [log10(HH8(:,n)); log10(HV8(:,n)); log10(VV8(:,n))];

        %Total u(:,n) for all datasets combined in the pixel
        up = u1(:,n);
        up(4:6) = u2(:,n);
        up(7:9) = u3(:,n);
        up(10:12) = u4(:,n);
        up(13:15) = u5(:,n);
        up(16:18) = u6(:,n);
        up(19:21) = u7(:,n);
        up(22:24) = u8(:,n);

        u(:,n) = up;

    %    dtotal(:,n) = (1/2)*transpose(u(:,n)-meanU)*inv(meanSigma)*(u(:,n)-meanU)+(1/2)*log(det(meanSigma));%
    y(:,n,k) = mvnpdf(u(:,n),meanU);
    % dtotal1(:,n) = (1/2)*transpose(u(:,n)-meanU);
    % dtotal2(:,n) = inv(meanSigma)*(u(:,n)-meanU);
    % dtotal3(:,n) = (1/2)*log(det(meanSigma));
    % dtotal(:,n) = dtotal1(:,n) + dtotal2(:,n) + dtotal3(:,n);
    end
    

end
    
    y1 = zeros([1,1048576]);
    y2 = zeros([1,1048576]);
    y3 = zeros([1,1048576]);
    y4 = zeros([1,1048576]);
    y5 = zeros([1,1048576]);
    y6 = zeros([1,1048576]);
    y7 = zeros([1,1048576]);
    y8 = zeros([1,1048576]);
    y9 = zeros([1,1048576]);
    y10 = zeros([1,1048576]);
    
    y1 = y(:,:,1);
    y2 = y(:,:,2);
    y3 = y(:,:,3);
    y4 = y(:,:,4);
    y5 = y(:,:,5);
    y6 = y(:,:,6);
    y7 = y(:,:,7);
    y8 = y(:,:,8);
    y9 = y(:,:,9);
    y10 = y(:,:,10);
%% Classifying the information
image = zeros([3,1,1048576]);

for n = 1:1048576
    
    %Getting the minimum distance
    maxY = y10(:,n);
    if maxY < y9(:,n)
        maxY = y9(:,n);
    end
    if maxY < y8(:,n)
        maxY = y8(:,n);
    end
    if maxY < y7(:,n)
        maxY = y7(:,n);
    end
    if maxY < y6(:,n)
        maxY = y6(:,n);
    end
    if maxY < y5(:,n)
        maxY = y5(:,n);
    end
    if maxY < y4(:,n)
        maxY = y4(:,n);
    end
    if maxY < y3(:,n)
        maxY = y3(:,n);
    end
    if maxY < y2(:,n)
        maxY = y2(:,n);
    end
    if maxY < y1(:,n)
        maxY = y1(:,n);
    end

    
    %Colour area allocation
    if maxY == y1(:,n)
        image(:,:,n) = [0 0 1]; %Blue-->Mask 1
    elseif maxY == y2(:,n) 
        image(:,:,n) = [245/255 222/255 179/255]; %Wheat 245,222,179-->Mask 2
    elseif maxY == y3(:,n)
        image(:,:,n) = [1 1 0]; %Yellow-->Mask 3
    elseif maxY == y4(:,n)
        image(:,:,n) = [1 0 1]; %Magenta-->Mask 4
    elseif maxY == y5(:,n)
        image(:,:,n) = [0 1 0]; %Green-->Mask 5
    elseif maxY == y6(:,n)
        image(:,:,n) = [1 215/255 0]; %Gold-->Mask 6
    elseif maxY == y7(:,n)
        image(:,:,n) = [0 100/255 0]; %Dark Green-->Mask 7
    elseif maxY == y8(:,n)
        image(:,:,n) = [0 1 1]; %Cyan-->Mask 8
    elseif maxY == y9(:,n)
        image(:,:,n) = [1 165/255 0]; %Orange-->Mask 9
    elseif maxY == y10(:,n)
        image(:,:,n) = [1 192/255 203/255]; %Pink 255,192,203-->Mask 10
    end
end
figure;
ImgClass = reshape(image,3,1024,1024);
I1 = permute(ImgClass,[2 3 1]);
I2 = flip(I1 ,1);
img = imrotate(I2,-90);
imshow(img);
colorbar;
% plane1 = ImgClass(1,:,:);
% plane2 = ImgClass(2,:,:);
% plane3 = ImgClass(3,:,:);
% rgbImage(:,:,1) = uint8(plane1);
% rgbImage(:,:,2) = uint8(plane2);
% rgbImage(:,:,3) = uint8(plane3);
% imshow(rgbImage);


