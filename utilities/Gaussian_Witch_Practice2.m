function [inside_vec] = Gaussian_Witch_Practice(A2,vertices)
% Use the vertices of the input polygon

% % Define the mean and covariance matrix for the Gaussian bivariate function
% mu = meanXY;
% sigma = covXY;
% % sigma(2,1)=0;
% % sigma(1,2)=0;
% % Create an anonymous function for the Gaussian bivariate function
% gauss_fun = @(x, y) gaussian_bivariate(x, y, mu, sigma);
inside_vec=zeros(size(A2,1),1);
for pp=1:1:size(A2,1)
point=A2(pp,:);
% Define the polygon mask function
polygon_mask = @(x, y) inpolygon(x, y, vertices(:,1), vertices(:,2));

% Create the integrand with the mask
% integrand = @(x, y) gauss_fun(x, y) .* polygon_mask(x, y);
% integrand = @(x, y) gauss_fun(x, y);
% Define the integration limits based on the polygon bounds

x_min = min(vertices(:,1));
x_max = max(vertices(:,1));
y_min = min(vertices(:,2));
y_max = max(vertices(:,2));

isInside = polygon_mask(point(1), point(2));

% Display result
if isInside
    inside_vec(pp)=1;
end

% if x_min<mu(1)-7*sqrt(sigma(1,1))
% x_min=mu(1)-7*sqrt(sigma(1,1));
% flag_x=1;
% else
% a=1;
% end
% 
% if x_max > mu(1)+7*sqrt(sigma(1,1)) || flag_x==1
% x_max = mu(1)+7*sqrt(sigma(1,1));
% else
% a=1;
% end
% 
% if y_min< mu(2)-7*sqrt(sigma(2,2))
% y_min= mu(2)-7*sqrt(sigma(2,2));
% flag_y=1;
% else
% a=1;
% end
% 
% if y_max > mu(2)+7*sqrt(sigma(2,2))  || flag_y==1
% y_max = mu(2)+7*sqrt(sigma(2,2));
% else
% a=1;
% end
% % Perform the numerical integration
% integral_value = integral2(integrand, x_min, x_max, y_min, y_max);
% Gaussian_error=1-integral_value;
% if Gaussian_error<0
% Gaussian_error=eps;
% end
% % Display the result
% fprintf('The integral error beyond the polygon is: %.4f\n', Gaussian_error);
% end
end
end
