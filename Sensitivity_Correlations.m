
%% This function performs a sensitivity analysis, by computing the correlation 
%% coefficients between the resoponse and the parameters and return a Pareto  plot.  

%% Author: Celine Scheidt


function s = Sensitivity_Correlations(real_params_val,Response,param_names)

%% Input Parameters
% - real_params_val: property values of the models (preferably defined by
%   ED). One column is one parameter.  Should be numerical values.
% - Response: vector of reponse of interest at one particular time
% - param_names: names of the parameters


%1. Compute the correlations of responses v.s paramater
 s = corr(real_params_val,Response);

%2. Create a Pareto plot
 [valsort,idxsort] = sort(abs(s));
 figure
 axes('FontSize',13);
 barh(abs(s(idxsort)))
 title('Sensitivity of parameters on Response','FontSize',15)
 set(gca,'YTickLabel',param_names(idxsort))
 

end
