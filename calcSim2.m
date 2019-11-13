function allcorr = calcSim2(spkmat,spkmat2,method,full)
% Sams script, but i think that it calculates either the dprime or the
% cosine or the mahalanobis distance... I think it would be a good idea to
% zscore this or normalize it before you do the calculations, so make sure
% that is done in the parent script
% the squareform for the p-dist function preserves the square format so
% that it can be indexed with a square matrix


if ~exist('full','var')
    full=2;
end
if isempty(full)
    full=2;
end


switch lower(method)
    
    case 'corr'
        % transposes spkmat, then takes the corr coef of every row to all
        % others
        
        %allcorr=corrcoef(spkmat,'rows','pairwise');
        allcorr=corr(spkmat',spkmat2','rows','pairwise','type','Spearman');
        % then takes the bottom triangle of that, because its a square
        % matrix
        %allcorr(tril(ones(size(allcorr)))==1) = nan;
    case 'mahalanobis'
        % takes the mahalanobis distance of the covariance matrix... it
        % must be in number of cells dimensions?
        allcorr= squareform(pdist2(spkmat,spkmat2,'mahalanobis'));
        %allcorr= squareform(pdist(spkmat(:,~all(C==0,1)),'mahalanobis'));
    case 'cosine'
        allcorr= pdist2(spkmat,spkmat2,method);
        allcorr=1-allcorr;
    case 'pearson'
        allcorr=corr(spkmat',spkmat2');
    case 'spearman'
        allcorr=corr(spkmat',spkmat2','rows','pairwise','type','Spearman');
        
    otherwise
        allcorr= squareform(pdist2(spkmat,spkmat2,method));
end
% kill lower half, or ust eye
if full>0
    % kill eye
    allcorr(logical(eye(length(allcorr))))=nan;
    if full>1
        % kill lower half
        allcorr(tril(ones(size(allcorr)))==1) = nan;
    end
    % kill eye
    
end
end
