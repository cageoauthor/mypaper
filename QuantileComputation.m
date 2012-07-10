%
% This function computes the quantiles of the values in response.

% Author: Celine Scheidt
% Date: April 2009


function Quantiles = QuantileComputation(Response,p,Clustering)

%% Input Parameters
%   - Response: matrix of the responses
%   - p: scalar or vector of cumulative probability values
%   - Clustering: results of the clustering (as given by function kernel_kmedoid). Optional 

%% Output Parameters 
%   - Quantiles: quantiles of the values in response.


Quantiles = zeros(size(Response,2),length(p));

if nargin == 2 % No clustering: each model as a weight of one
    for i = 1:size(Response,2)
        Quantiles(i,:)=quantile(Response(:,i),p);
    end
else % Each model is weighted by the number of point in the centroid (weights)
    WeigthedResponse = [];
    for i = 1:size(Clustering.medoids,2)
        WeigthedResponse=vertcat(WeigthedResponse,repmat(Response(Clustering.medoids(i),:),Clustering.weights(i),1));
    end
    for i = 1:size(Response,2)
        Quantiles(i,:)=quantile(WeigthedResponse(:,i),p);
    end
    
    
end
    
end