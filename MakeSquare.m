function matout = MakeSquare(matin,cutwhere)
% no longer pads, but instead cuts, this makes your matrix a square, is
% especially useful if you have coordinates that are tracking out of your
% maze i
% it cuts equally on both sides, one good addition may be to add an option
% to cut left, top, right, or bottom.
% cutwhere: if cutwhere =0, it cuts evenly top, bottom, left, right
%  if cutwhere is 1, it cuts right or bottom
% if cutwhere is 2, it cuts left, or top


if ~exist('cutwhere','var')
    cutwhere=0;
end


a = diff(size(matin));
if cutwhere==0
    if a<0 % if there are more rows than cols
        
        matout = matin(1+floor(-a/2):end-ceil(-a/2), :);
        
    elseif a>0 % if there are more cols than rows
        
        matout = matin(:, 1+floor(a/2):end-ceil(a/2));
        
    else
        matout = matin;
        
    end
elseif cutwhere==1
    if a<0 % if there are more rows than cols
        
        matout = matin(1:floor(size(matin)), :);
        
    elseif a>0 % if there are more cols than rows
        
        matout = matin(:, 1:floor(size(matin)));
        
    else
        matout = matin;
        
    end
elseif cutwhere==2
    if a<0 % if there are more rows than cols
        
        matout = matin(abs(a):end, :);
        
    elseif a>0 % if there are more cols than rows
        
        matout = matin(:, abs(a):end);
        
    else
        matout = matin;
        
    end
end