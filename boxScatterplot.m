function boxScatterplot(x,grps,varargin)
%boxScatterplot(x,grps,varargin)
%
%   I usually like to plot individual data points on top of boxplots so
%   this function does that by separating your x vector into grps and then
%   jittering the data points to overlay. 
%
%   INPUTS
%       x: vector of data points to plot.
%
%       grps: grouping vector the same size as x. Each value indicates
%       which group each value of x belongs to. 
%
%       NAME,VALUE arguments:
%           xLabels: cell array of strings specifying the group names, must
%           be in the same order as grps (i.e., the min of grps must be the
%           first string in xLabels and the max of grps must be the last. 
%           
%           yLabel: string specifying y axis label. 
%
%           boxColor: string or RGB value specifying what color you want
%           the boxplot to be.
%
%           circleSize: scalar or vector specifying scatter plot circle
%           sizes. 
%
%           circleColors: string or RGB values specifying what color you
%           want all or individual circles in the scatter plot to be.
%       
%           transparency: scalar specifying how transparent you want the
%           circles to be. 1 is opaque, 0 is completely transparent. 
%
%           sf: spread factor, scalar specifying how wide you want the
%           individual data points to be jittered. 
%
%           position: vector specifying the position of the figure. 
%
%           plotBox: whether or not to make boxplot.
%

%% Set up.
    p = inputParser; 
    p.addRequired('x',@(x) isnumeric(x));
    p.addRequired('grps',@(x) isnumeric(x));
    % will forgot to match the labels with the number of groups... asshole

    p.addParameter('xLabels',{'Group 1','Group 2'},@(x) iscell(x));
    p.addParameter('yLabel','Metric',@(x) ischar(x)); 
    p.addParameter('boxColor','k',@(x) ischar(x) || isnumeric(x));
    p.addParameter('circleSize',20,@(x) isnumeric(x)); 
    p.addParameter('circleColors',[.5 .5 .5],@(x) ischar(x) || isnumeric(x));
    p.addParameter('transparency',.3,@(x) isscalar(x)); 
    p.addParameter('sf',.05,@(x) isscalar(x));
    p.addParameter('position',[520 350 300 450]); 
    p.addParameter('plotBox',true,@(x) islogical(x));
    
    p.parse(x,grps,varargin{:});
    xLabels = p.Results.xLabels; 
    yLabel = p.Results.yLabel; 
    boxColor = p.Results.boxColor;
    circleSize = p.Results.circleSize;
    transparency = p.Results.transparency;
    circleColors = p.Results.circleColors;
    sf = p.Results.sf; 
    position = p.Results.position;
    plotBox = p.Results.plotBox;
    
    % fix xlabels
    
    if length(unique(grps))>length(xLabels)
        for i=1:length(unique(grps))
            xLabels{i}=sprintf('group %d',i);
        end
    end
    
    % put values and groups into a design matrix
    try
        design(:,1)=x(:); design(:,2)=grps(:);
        design=sortrows(design,2);
    catch
        fprintf('x and grps have to be the same size');
    end
    
    %Number of values. 
    n = size(design,1);
%% Make figure. 
    %Group numbers. 
    grpNums = unique(design(:,2))';
    jitters = zeros(n,1); 
    
    %Create jitter.
    c = 1; 
    for g = grpNums
        %Number of elements corresponding to that group. 
        nInGrp = sum(design(:,2)==g);
        
        %Jitter!
        jitters(c:c+nInGrp-1) = g - (sf*randn(nInGrp,1));
        
        %Step. 
        c = c+nInGrp; 
    end
    
    %Figure here. 
    if length(position)==4
        figure('Position',position); 
    end
    hold on;

    if plotBox
        boxplot(design(:,1),design(:,2),'color',boxColor,'symbol','k',...
            'positions',unique(design(:,2)));
        boxProps = get(gca,'Children');
        [boxProps(1).Children.LineWidth] = deal(2);
    end
    
    scat = scatter(jitters,design(:,1),circleSize,circleColors,'filled');
    alpha(scat,transparency);
    
    ylabel(yLabel);
    set(gca,'tickdir','out','XTickLabels',xLabels);
   
    
    
end