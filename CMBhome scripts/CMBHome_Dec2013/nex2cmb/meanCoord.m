function [t,x,y] = meanCoord(c)
    t = c(:,1);
    if(size(c,2)==4); c = c(:,1:3); end % Remove 'motion'
    
    if(size(c,2)>=3); c((c(:,2)==0 & c(:,3)==0),2:3) = NaN; end
    if(size(c,2)>=5); c((c(:,4)==0 & c(:,5)==0),4:5) = NaN; end
    if(size(c,2)>=7); c((c(:,6)==0 & c(:,7)==0),6:7) = NaN; end
    
    x = nanmean(c(:,2:2:end),2);
    y = nanmean(c(:,3:2:end),2);

    x(isnan(x)) = 0;
    y(isnan(y)) = 0;
    
    if(nargout == 1); t = [t x y]; end
end