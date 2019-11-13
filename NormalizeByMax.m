function [normalized] = NormalizeByMax(matrix,dim)
% [normalized] = NormalizeByMax(matrix,dim)
% normalizes each column or row by the max of that col/row
% dim 1 is normalize across first dim, 2 is normalize across 2nd (row)

if ~exist('dim','var')
    dim=1;
end


switch dim
    case 1 % normalize each column to 1
        for i=1:size(matrix,2)
            matrix(:,i)=(matrix(:,i)./max(matrix(:,i)));
        end
    case 2 % each row by 1
        for i=1:size(matrix,1)
            matrix(i,:)=(matrix(i,:)./max(matrix(i,:)));
        end
    case 3 % each depth by 1
        try
            for i=1:size(matrix,1)
                for j=1:size(matrix(2))
                    matrix(i,j,:)=(matrix(i,j,:)./max(matrix(i,j,:)));
                end
            end
        catch
            fprintf('Your matrix is only two dimensions!');
        end
end
normalized=matrix;
end

