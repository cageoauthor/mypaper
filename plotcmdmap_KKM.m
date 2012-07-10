function plotcmdmap_KKM(Xd,KKM)

%% Function which maps the points defined in Xd and (optional) color the
%% points according to clustering results obtained by kmedoid.m
% Celine Scheidt

% Input Parameters:
% - Xd: location of the points to be plotted (from MDS). Rows of Xd are the coordinates of n points in p-dimensional space
% - KKM: results from the clustering in the format givent by kmedoid. This is optional.  
    % If given, point are colored by clusters; medoids are represented by squares 

if nargin ==2
    Dmax_ = max(KKM.T);
    Dmin_ = min(KKM.T);
    Dmax = Dmax_ + 0.02 * (Dmax_ - Dmin_) + 0.0001;
    Dmin = Dmin_ - 0.02 * (Dmax_ - Dmin_) - 0.0001;
    Cs = [1 0 1;1 0 0;1 1 0;0 1 0;0 1 1;0 0 1];
    Ds = linspace(Dmin,Dmax,size(Cs,1));
    C = interp1(Ds', Cs, unique(KKM.T), 'linear');

    figure
    for i=1:length(KKM.medoids)
        plot(Xd(find(KKM.T==i),1), Xd(find(KKM.T==i),2), 'o', 'Color',C(i,:), ...
            'MarkerSize', 6, 'LineWidth', 1,'MarkerFaceColor',C(i,:), ...
            'MarkerEdgeColor','k');

        hold on
        plot(Xd(KKM.medoids(i),1), Xd(KKM.medoids(i),2), 'bs', ...
            'MarkerSize', 10, 'LineWidth', 2,'MarkerFaceColor',C(i,:), ...
            'MarkerEdgeColor','k');
    end

else
    C = zeros(size(Xd,1),3);
    for i = 1:size(Xd,1)
        C(i,:)=[1 0 0];
    end
    
    figure
    plot(Xd(:,1), Xd(:,2), 'o', 'Color',C(i,:), ...
            'MarkerSize', 6, 'LineWidth', 1,'MarkerFaceColor',C(i,:), ...
            'MarkerEdgeColor','k');
end


hold off
grid on

set(gca,'LineWidth',3)
set(gca,'XTickLabel',{''})
set(gca,'YTickLabel',{''})

end
