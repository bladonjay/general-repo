function [smoothmat] = GaussNdim(kernel_size,std)
% function [smoothmat] = GaussNdim(kernel_size,std)
% generates a smoothing kernel in n dimensions
% INPUTS:
% kernel size in n dimensions
% std: the bin with of the standard deviation
% 10 december 2009 andrew 
% this sequel uses bens nanconvolve script
% JHB 9-4-14

% force standard dev to be bigger than 0
if nargin<2
    std=1;
end
if std == 0, return; end
% force even dimension sizes
if any(rem(kernel_size,2))
    warning('some dimension inputs sizes were odd, adding to make them even');
    kernel_size=kernel_size+rem(kernel_size,2);
end
    


if numel(kernel_size)==1
    grid=meshgrid(-kernel_size/2:kernelsize/2);
    kernel = pdf('Normal', grid, 0, std);
    smoothmat = kernel./sum(sum(kernel));
elseif numel(kernel_size)==2
    
    [Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);

    Rgrid=sqrt((Xgrid.^2+Ygrid.^2));

    kernel = pdf('Normal', Rgrid, 0, std);

    smoothmat = kernel./sum(sum(kernel));
    
elseif numel(kernel_size)==3
    [Xgrid,Ygrid,Zgrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2, -kernel_size(3)/2:kernel_size(3)/2);

    Rgrid=sqrt((Xgrid.^2+Ygrid.^2+Zgrid.^2));

    % now linearize
    
    kernel = pdf('Normal', Rgrid(:), 0, std);
    % and reshape

    smoothvector = kernel./sum(kernel);
    smoothmat=reshape(smoothvector,[kernel_size(1)+1,kernel_size(2)+1,kernel_size(3)+1]);
else
    error('cant do more than three dimensions');
end
   
end

