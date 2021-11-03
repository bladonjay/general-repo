function mat = SmoothMat2(mat, kernel_size, std, nanout)
%
% Smooths matrix by convolving with 2d gaussian of size
% kernel_size=[bins_x bins_y] and standard deviation 'std'
%
% if std==0, just returns mat
%
% 10 december 2009 andrew
% this sequel uses bens nanconv script
% JHB 9-4-14

if ~exist('nanout','var')
    nanout=true;
end

if nargin<3
    std=1;
end

if std == 0, return; end

% its smart to force both dimensions to be odd.
[Xgrid,Ygrid]=meshgrid(floor(-kernel_size(1)/2): ceil(kernel_size(1)/2),...
    floor(-kernel_size(2)/2):ceil(kernel_size(2)/2));

if length(std)==1
    
    Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
    
    kernel = pdf('Normal', Rgrid, 0, std);
    
    kernel = kernel./sum(sum(kernel));
    
elseif length(std)==2

    Rgrid=sqrt(((Xgrid./std(1)).^2+(Ygrid./std(2)).^2));

    kernel = pdf('Normal', Rgrid, 0, 1);
end
mat = nanconv(mat, kernel, 'same','edge','nanout',nanout);
end

