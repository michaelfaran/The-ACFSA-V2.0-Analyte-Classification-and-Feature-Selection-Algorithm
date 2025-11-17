function [avgErr, classErr, cutpoints, zoneLabels] = qda1d_gaussians(mu_vec, var_vec)
% QDA in 1D for K Gaussians (equal priors).
% Inputs:
%   mu_vec  : 1xK or Kx1 vector of means (can be any order)
%   var_vec : 1xK or Kx1 vector of variances (>0)
% Outputs:
%   avgErr     : scalar average misclassification error (equal priors)
%   classErr   : Kx1 per-class misclassification probabilities
%   cutpoints  : sorted unique boundaries including -Inf and +Inf
%   zoneLabels : class label for each interval (cutpoints(i), cutpoints(i+1))
%
% Notes:
% - True 1D QDA may produce two cutpoints between a pair if variances differ.
% - We merge all pairwise cutpoints and classify intervals by argmax g_k(x).

    mu   = mu_vec(:).';    % make row
    sig2 = var_vec(:).';   % make row
    K = numel(mu);
    if K < 2, error('Need at least 2 classes.'); end
    if K >= 10
        warning('Function intended for K<10; continuing with K=%d.', K);
    end
    if any(sig2 <= 0)
        error('All variances must be positive.');
    end

    % ---- collect all pairwise QDA boundaries ----
    roots_all = [];
    for i = 1:K-1
        for j = i+1:K
            r = qda_pair_boundaries_1d(mu(i), sig2(i), mu(j), sig2(j));
            if ~isempty(r)
                roots_all = [roots_all; r(:)]; %#ok<AGROW>
            end
        end
    end
    roots_all = roots_all(isfinite(roots_all) & isreal(roots_all));
    if ~isempty(roots_all)
        roots_all = unique_tol(roots_all, 1e-12);
    end

    % ---- final cutpoints and interval labels ----
    cutpoints = [-Inf; roots_all(:); Inf];
    M = numel(cutpoints)-1;
    zoneLabels = zeros(M,1);

    for m = 1:M
        a = cutpoints(m); b = cutpoints(m+1);
        if isinf(a) && ~isinf(b)
            xmid = b - min(1, abs(b)+1);
        elseif ~isinf(a) && isinf(b)
            xmid = a + min(1, abs(a)+1);
        elseif isinf(a) && isinf(b)
            xmid = 0;
        else
            xmid = 0.5*(a+b);
        end
        g = qda_discriminants_1d(xmid, mu, sig2);   % 1xK
        [~, zoneLabels(m)] = max(g);
    end

    % ---- per-class correct probs & errors (analytic integration) ----
    correctPerClass = zeros(K,1);
    for k = 1:K
        pk = 0;
        for m = 1:M
            if zoneLabels(m) == k
                a = cutpoints(m); b = cutpoints(m+1);
                pk = pk + normal_interval_prob(a, b, mu(k), sqrt(sig2(k)));
            end
        end
        correctPerClass(k) = pk;
    end
    classErr = 1 - correctPerClass;                % Kx1
    avgErr   = 1 - mean(correctPerClass);          % equal priors ⇒ mean

end

% ---------- helpers ----------
function r = qda_pair_boundaries_1d(m1, s2_1, m2, s2_2)
% Solve for g1(x) = g2(x) in 1D QDA (equal priors):
%  (x-m1)^2/(2 s2_1) + 0.5*log s2_1 = (x-m2)^2/(2 s2_2) + 0.5*log s2_2
% => A x^2 + B x + C = 0
    A = (1/s2_1) - (1/s2_2);
    B = -2*(m1/s2_1 - m2/s2_2);
    C = (m1^2/s2_1) - (m2^2/s2_2) + log(s2_1/s2_2);
    tol = 1e-14;
    if abs(A) < tol
        % nearly equal variances → linear
        if abs(B) < tol
            r = [];  % degenerate: indistinguishable
        else
            r = -C / B;
        end
    else
        D = B.^2 - 4*A*C;
        if D < 0
            r = [];  % no real intersections
        else
            sqrtD = sqrt(max(D,0));
            r = [(-B - sqrtD)/(2*A), (-B + sqrtD)/(2*A)];
        end
    end
end

function g = qda_discriminants_1d(x, mu, sig2)
% g_k(x) up to common additive constant (equal priors)
% g_k(x) = -0.5*((x - mu_k)^2/sig2_k + log(sig2_k))
    g = -0.5*(((x - mu).^2)./sig2 + log(sig2));
end

function p = normal_interval_prob(a, b, m, s)
% P(a < X < b) for X ~ N(m, s^2)
    if isinf(a) && a < 0, Fa = 0; else, Fa = normcdf((a - m)/s); end
    if isinf(b) && b > 0, Fb = 1; else, Fb = normcdf((b - m)/s); end
    p = max(0, min(1, Fb - Fa));
end

function u = unique_tol(v, tol)
% Unique with tolerance (ascending)
    v = sort(v(:));
    if isempty(v), u = v; return; end
    u = v(1);
    for i = 2:numel(v)
        if abs(v(i) - u(end)) > tol
            u(end+1,1) = v(i); %#ok<AGROW>
        end
    end
end
