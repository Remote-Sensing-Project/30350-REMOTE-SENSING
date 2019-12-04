%% Confusion Matrix calculation
 classID = zeros([1,1048576]);

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
        classID(:,n) = 1;
    elseif maxY == y2(:,n) 
        classID(:,n) = 2;
    elseif maxY == y3(:,n)
        classID(:,n) = 3;
    elseif maxY == y4(:,n)
        classID(:,n) = 4;
    elseif maxY == y5(:,n)
        classID(:,n) = 5;
    elseif maxY == y6(:,n)
        classID(:,n) = 6;
    elseif maxY == y7(:,n)
        classID(:,n) = 7;
    elseif maxY == y8(:,n)
        classID(:,n) = 8;
    elseif maxY == y9(:,n)
        classID(:,n) = 9;
    elseif maxY == y10(:,n)
        classID(:,n) = 10;
    end
end


Nclasses = 10;
confusion = zeros([10,10]);
for k = 1:Nclasses
    q = find(masks(:,:,k) == 1);
    classIDmask = classID(:,q);
        for pix = 1:length(q)
            if classIDmask(:,pix)== 1
                confusion(k,1) = confusion(k,1) + 1;
            end
            if classIDmask(:,pix)== 2
                confusion(k,2) = confusion(k,2) + 1;
            end
            if classIDmask(:,pix)== 3
                confusion(k,3) = confusion(k,3) + 1;
            end
            if classIDmask(:,pix)== 4
                confusion(k,4) = confusion(k,4) + 1;
            end
            if classIDmask(:,pix)== 5
                confusion(k,5) = confusion(k,5) + 1;
            end
            if classIDmask(:,pix)== 6
                confusion(k,6) = confusion(k,6) + 1;
            end
            if classIDmask(:,pix)== 7
                confusion(k,7) = confusion(k,7) + 1;
            end
            if classIDmask(:,pix)== 8
                confusion(k,8) = confusion(k,8) + 1;
            end
            if classIDmask(:,pix)== 9
                confusion(k,9) = confusion(k,9) + 1;
            end
            if classIDmask(:,pix)== 10
                confusion(k,10) = confusion(k,10) + 1;
            end
        end
end

%% Average of the diagonal
mDiagTr = trace(confusion)/10;
diagonal = diag(confusion);
q = zeros([1 10]);
acc = zeros([1 10]);
for k = 1:Nclasses
    q = find(masks(:,:,k) == 1);
    acc(k) = diagonal(k)/length(q);
end
accuracy = mean(acc);
%%
cm = confusionchart(confusion);
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';