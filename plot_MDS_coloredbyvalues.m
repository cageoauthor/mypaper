
%% This function colors points in MDS according to a particular value 
 
% Author: Celine Scheidt
% Date: Nov. 2010

function plot_MDS_coloredbyvalues(Xd,values)

%% Input Parameters:
%    - Xd: Location of the points in the MDS map
%    - values:  Vector of values at each point



Rmax_ = max(values);
Rmin_ = min(values);
Rmax = Rmax_ + 0.02 * (Rmax_ - Rmin_) + 0.0001;
Rmin = Rmin_ - 0.02 * (Rmax_ - Rmin_) - 0.0001;

Cs = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1];

Rs = linspace(Rmin,Rmax,size(Cs,1));
C = interp1(Rs', Cs, values, 'linear');

figure
for i=1:size(Xd,1)
    plot(Xd(i,1), Xd(i,2), 'o', 'Color',C(i,:), ...
        'MarkerSize', 6,'LineWidth',1,'MarkerFaceColor',C(i,:),...
        'MarkerEdgeColor','k');
    hold on
end
grid on;
set(gca,'LineWidth',3)
set(gca,'XTickLabel',{''})
set(gca,'YTickLabel',{''})


end
