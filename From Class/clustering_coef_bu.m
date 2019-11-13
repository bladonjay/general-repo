function C=clustering_coef_bu(G)
%CLUSTERING_COEF_BU     Clustering coefficient
%
%   C = clustering_coef_bu(A);
%
%   The clustering coefficient is the fraction of triangles around a node
%   (equiv. the fraction of node’s neighbors that are neighbors of each other).
%
%   Input:      G,      binary undirected connection matrix
%
%   Output:     C,      clustering coefficient vector
%
%   Reference: Watts and Strogatz (1998) Nature 393:440-442.
%
%
%   Mika Rubinov, UNSW, 2007-2010

% number of nodes
n=length(G);
% our output
C=zeros(n,1);

% for each row
for u=1:n
    % find where that node has neighbors
    V=find(G(u,:));
    % for each element, get the 1st neighbors
    k=length(V);
    if k>=2;                %degree must be at least 2
        % find how many of its 1st neighbors connect to eachother
        S=G(V,V);
        % k*k-1 is total number of neighbors
        C(u)=sum(S(:))/(k^2-k);
        
    end
end
% end up with the proportion of neighbors that are connected