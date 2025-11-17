function [handle_ellipse,meanXY,covXY] = plot_ellipse_inflate(x,y,nn,STD_buff,Inflate)    
% handle_ellipse - simply plots an ellipse arround scattered points with 95% confidence interval
% The code based on Ohad Gal fit_ellipse (2003)

%
% Format:   handle_ellipse  = plot_ellipse( x,y )    
%
% Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
%
% Output:   handle_plot - handle of the plot object
%
% Note:     This is only simple function without amenities.
%edited by Michael to output the mu and sigma

%inflate factor due to limited measrument
n=nn;
if STD_buff==0 && Inflate==1
  if n <= 1
        error('Need at least 2 samples for a variance estimate.');
    end

    nu = n;  % degrees of freedom

    if n > 3
        faktor = ((nu-1) / (nu - 2)) * ((n+1)/n) ; %factor for the varaince
    else
        % Use the posterior median instead
        % Quantile function for Inv-χ²: sigma^2 = nu*s2 / χ²_{nu,1-p}
        % with p = 0.5 (median)
        q = chi2inv(0.5, nu);         % median of χ²_{nu}
        faktor = (nu *((n+1)/n)) / q;
    end
elseif STD_buff>0 && Inflate==0
faktor=STD_buff^2;
else
   if n <= 1
        error('Need at least 2 samples for a variance estimate.');
    end

    nu = n;  % degrees of freedom

    if n > 2
        faktor = ((nu-1) / (nu - 2)) * ((n+1)/n) ; %factor for the varaince
    else
        % Use the posterior median instead
        % Quantile function for Inv-χ²: sigma^2 = nu*s2 / χ²_{nu,1-p}
        % with p = 0.5 (median)
        q = chi2inv(0.5, nu);         % median of χ²_{nu}
        faktor = (nu *((n+1)/n)) / q;
    end
    faktor=faktor.*STD_buff^2;
end
% Center of the ellipse
    X0=mean(x);
    Y0=mean(y);
    covXY=cov(x,y).*faktor;
    meanXY=[X0 Y0];
% Lenght of the axes of the ellipse
    x=x-mean(x);
    y=y-mean(y);
    l= PCA([ x,y]);
    orientation_rad =-atan(l(2,1)/l(1,1));
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    xy=[x,y]*R;
    a=2*std(xy(:,1)).*faktor;
    b=2*std(xy(:,2)).*faktor;

% Matrix of rotation
    R = [ cos_phi -sin_phi; sin_phi cos_phi ];

% Drawing the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     =  a*cos( theta_r );
    ellipse_y_r     = b*sin( theta_r );
    rotated_ellipse =  [ellipse_x_r;ellipse_y_r]'*R;
    rotated_ellipse(:,1)=rotated_ellipse(:,1) + X0;
    rotated_ellipse(:,2)=rotated_ellipse(:,2) + Y0;
    handle_ellipse =plot( rotated_ellipse(:,1),rotated_ellipse(:,2),'r' );
end


function loadings=PCA(X)
    D=bsxfun(@minus, X, mean(X, 1));
    CVD=D'*D/size(X,1);
    [~,I]=sort(eig(CVD),'descend');
    [H,V]=eig(CVD);
    loadings=H(:,I);
end