function [Gaussian_error] = Gaussian_Witch_Practice(meanXY,covXY,vertices)
% Use the vertices of the input polygon

% Define the mean and covariance matrix for the Gaussian bivariate function
mu = meanXY;
sigma = covXY;
% sigma(2,1)=0;
% sigma(1,2)=0;
% Create an anonymous function for the Gaussian bivariate function
gauss_fun = @(x, y) gaussian_bivariate(x, y, mu, sigma);

% Define the polygon mask function
polygon_mask = @(x, y) inpolygon(x, y, vertices(:,1), vertices(:,2));

% Create the integrand with the mask
integrand = @(x, y) gauss_fun(x, y) .* polygon_mask(x, y);
% integrand = @(x, y) gauss_fun(x, y);
% Define the integration limits based on the polygon bounds

x_min = min(vertices(:,1));
x_max = max(vertices(:,1));
y_min = min(vertices(:,2));
y_max = max(vertices(:,2));
flag_x=0;
flag_y=0;
if x_min<mu(1)-7*sqrt(sigma(1,1))
x_min=mu(1)-7*sqrt(sigma(1,1));
flag_x=1;
else
a=1;
end

if x_max > mu(1)+7*sqrt(sigma(1,1)) || flag_x==1
x_max = mu(1)+7*sqrt(sigma(1,1));
else
a=1;
end

if y_min< mu(2)-7*sqrt(sigma(2,2))
y_min= mu(2)-7*sqrt(sigma(2,2));
flag_y=1;
else
a=1;
end

if y_max > mu(2)+7*sqrt(sigma(2,2))  || flag_y==1
y_max = mu(2)+7*sqrt(sigma(2,2));
else
a=1;
end
% Perform the numerical integration
integral_value = integral2(integrand, x_min, x_max, y_min, y_max);
Gaussian_error=1-integral_value;
if Gaussian_error<0
Gaussian_error=eps;
end
% Display the result
fprintf('The integral error beyond the polygon is: %.4f\n', Gaussian_error);
end


function z = gaussian_bivariate(x, y, mu, sigma)
    % x, y are coordinates
    % mu is the mean vector [mu_x, mu_y]
    % sigma is the covariance matrix [sigma_xx, sigma_xy; sigma_yx, sigma_yy]
    mu_x = mu(1);
    mu_y = mu(2);
    sigma_xx = sigma(1,1);
    sigma_xy = sigma(1,2);
    sigma_y = sigma(2,2);
    rho=sqrt(sigma_xy.^2./(sigma_xx*sigma_y));
    % Gaussian bivariate formula
    z = exp(-0.5*(1/(1-rho^2)) * ((x - mu_x).^2 / sigma_xx + (y - mu_y).^2 / sigma_y - 2*sigma_xy*(x - mu_x).*(y - mu_y) / (sigma_xx * sigma_y))) ...
        / (2 * pi * sqrt(det(sigma)));
end
