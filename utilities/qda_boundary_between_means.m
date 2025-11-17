function x_star = qda_boundary_between_means(m1, m2, sigma1, sigma2)
%QDA_BOUNDARY_BETWEEN_MEANS 
%   Returns the decision boundary between two 1D Gaussians N(m1,sigma1^2) and N(m2,sigma2^2)
%   with equal priors. Keeps only the root that lies between the two means.
%
%   Usage:
%   x_star = qda_boundary_between_means(m1, m2, sigma1, sigma2)

    s1sq = sigma1^2;
    s2sq = sigma2^2;
    delta_m = m2 - m1;
    delta_s2 = s2sq - s1sq;

    % Equal variances -> midpoint
    if abs(delta_s2) < 1e-12
        x_star = (m1 + m2) / 2;
        return;
    end

    % Term under the square root
    term = (delta_m)^2 + 2*delta_s2*log(sigma2/sigma1);
    if term < 0
        x_star = NaN; % no real solutions
        return;
    end

    sqrt_term = sigma1 * sigma2 * sqrt(term);

    numerator_plus  = s2sq*m1 - s1sq*m2 + sqrt_term;
    numerator_minus = s2sq*m1 - s1sq*m2 - sqrt_term;
    denom = delta_s2;

    x_plus  = numerator_plus / denom;
    x_minus = numerator_minus / denom;

    % Keep only the root that lies between m1 and m2
    lower = min(m1, m2);
    upper = max(m1, m2);
    roots_between = [x_plus, x_minus];
    mask = (roots_between >= lower) & (roots_between <= upper);

    if any(mask)
        x_star = roots_between(mask);
        % If both somehow inside, pick the one closer to midpoint
        if numel(x_star) > 1
            x_star = x_star( abs(x_star - (m1+m2)/2) == min(abs(x_star - (m1+m2)/2)) );
        end
    else
        x_star = NaN; % none inside
    end
end
