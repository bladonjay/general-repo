function mat = SmoothMat2(mat, kernel_size, std)
%
% Smooths matrix by convolving with 2d gaussian of size
% kernel_size=[bins_x bins_y] and standard deviation 'std'
%
% if std==0, just returns mat
%
% 10 december 2009 andrew
% this sequel uses bens nanconv script
% JHB 9-4-14

if nargin<3
    std=1;
end

if std == 0, return; end

[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);

Rgrid=sqrt((Xgrid.^2+Ygrid.^2));

kernel = pdf('Normal', Rgrid, 0, std);

kernel = kernel./sum(sum(kernel));

mat = nanconv(mat, kernel, 'same','edge','nanout',true);
end

