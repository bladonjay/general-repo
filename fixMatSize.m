function [matrix] = fixMatSize(matrix,sizes)
% this adds nans or removes extra values to reshape the 1 or 2d matrix to
% shape sizes
for i=1:length(sizes)
    if i==1
        if size(matrix,i)>sizes(i)
            matrix=matrix(1:sizes(i),:);
        elseif size(matrix,i)<sizes(i)
            matrix(size(matrix,1)+1:sizes(i),:)=nan;
        end
    elseif i==2
        if size(matrix,i)>sizes(i)
            matrix=matrix(:,1:sizes(1));
        elseif size(matrix,i)<sizes(i)
            matrix(:,size(matrix,1)+1:sizes(i))=nan;
        end
    end
end