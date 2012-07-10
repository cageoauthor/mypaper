
%% This function applies a kernel k-medoid algorithm (clustering) 
%% in the Feature Space defined by the kernel matrix K 
 
% Author: Celine Scheidt
% Date: April 2009


function Clustering = kernel_kmedoid(Xd,nbclusters,K)

%% Input Parameters
%   - Xd: matrix (N_models x Dim_MDS) of the coordinates of initial models
%   - nbclusters: number of clusters to be constructed
%   - K: Kernel matrix 

%% Output Parameters 
%   - Clustering: Results of the clustering, which contains:
%        - T: is a vector of length N_models, which contains the clusters each realization belongs to
%        - medoids: is a vector of length nbcluster containing the index
%        of the medoids
%        - weights: is a vector of ength nbcluster containing the number
%         of models in each cluster.
%


maxIterations = 50;
npoints = size(Xd,1);
%rand('seed',734636);
rand('seed',734600);

%% 1. Random selection of the inital nbclusters medoids.

initMedoids = randperm(npoints);
initMedoids = initMedoids(1:nbclusters);


%% 2. Associate each points to the closest medoid

% 2.1 Compute distance between each point and the selected initial medoids. Note that the distance is in the 
% feature space, so it can be computed using the kernel matrix only) 

dist_points_medoids = zeros(npoints,nbclusters);

for i = 1:npoints
    for j = 1:nbclusters
        dist_points_medoids(i,j) = K(i,i) - 2*K(i,initMedoids(j)) + K(initMedoids(j),initMedoids(j));
    end
end

% 2.2 Minimization: T contains the value of the cluster each model belongs to
[B,T] = min(dist_points_medoids,[],2);


%% 3.  Update the medoid of each cluster and re-assign each point to a cluster

currentMedoids_prev_iter = ones(1,nbclusters);
currentMedoids = initMedoids;
nbIter = 0;


while (~all(currentMedoids == currentMedoids_prev_iter) && nbIter < maxIterations)  % while cluster configuration is changing and maxIteration not reached
    currentMedoids_prev_iter = currentMedoids;
    
    % For each cluster
    for i = 1:nbclusters
        pts_in_cluster = find(T == i);
         % Compute distance between each point of the cluster and its medoid.
        dist_within_cluster = zeros(size(pts_in_cluster,2),size(pts_in_cluster,2));
        for j = 1:size(pts_in_cluster,1)
                for k = 1:size(pts_in_cluster,1)
                    dist_within_cluster(j,k) = K(pts_in_cluster(j),pts_in_cluster(j)) - 2*K(pts_in_cluster(j),pts_in_cluster(k)) + K(pts_in_cluster(k),pts_in_cluster(k));
                end
        end
        % minimize the distance and select the new medoid
        [dclust idx_min] = min(mean(dist_within_cluster));
        currentMedoids(i) = pts_in_cluster(idx_min);
    end

    % New medoids are defined, compute distance between the points and medoids
    for i = 1:npoints
        for j = 1:nbclusters
            dist_points_medoids(i,j) = K(i,i) - 2*K(i,currentMedoids(j)) + K(currentMedoids(j),currentMedoids(j));
        end
    end
        
    %Associate each points to the closest medoid
    [B,T] = min(dist_points_medoids,[],2);
    
    nbIter = nbIter +1;    
        
end

%% Once the medoids are defined, store the outputs

weights = zeros(nbclusters,1);
for i = 1:nbclusters
    weights(i) = sum(T == i);
end
    

Clustering.T = T;  
Clustering.medoids = currentMedoids;  
Clustering.weights = weights;

end