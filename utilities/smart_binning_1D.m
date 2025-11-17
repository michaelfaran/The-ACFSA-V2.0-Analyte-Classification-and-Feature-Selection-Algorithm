function [Binned, zoneInfo, mu_sorted, sigma_sorted] = smart_binning_1D(DAT_in, n_measure, all_analyte_length, Decision_boundaries)
% SMART_BINNING_1D  Per-sensor "smart" binning using 1D Gaussian decision zones.
%   [Binned, zoneInfo, mu_sorted, sigma_sorted] = smart_binning_1D(DAT_in, n_measure, K, Decision_boundaries)
%
% Inputs
%   DAT_in               : either a numeric matrix [N x P] or a .mat filename containing variable DAT
%   n_measure            : number of measurements per class (block size)
%   all_analyte_length   : number of classes K (blocks)
%   Decision_boundaries  : 1 -> Voronoi (nearest-mean midpoints), 0 -> QDA (1D Gaussian boundaries)
%
% Assumptions
%   - Samples are ordered by class in DAT: rows 1:n_measure -> class 1, next n_measure -> class 2, ..., class K.
%
% Outputs
%   Binned       : [N x P] integer labels in {1..K}, predicted by per-sensor decision zones
%   zoneInfo     : struct array (1 x P) with fields:
%                  .order (1xK original-class indices sorted by mean),
%                  .boundaries (1x(K-1) decision thresholds, ascending),
%                  .mu (1xK sorted means), .sigma (1xK sorted stds)
%   mu_sorted    : [K x P] sorted means per sensor (ascending order)
%   sigma_sorted : [K x P] stds reordered to match mu_sorted
%
% Example
%   [Binned, Z] = smart_binning_1D('myDAT.mat', 10, 5, 0);  % QDA zones

    % ---- Load DAT ----
    if ischar(DAT_in) || isstring(DAT_in)
        S = load(DAT_in);
        if isfield(S, 'DAT')
            DAT = S.DAT;
        else
            % take the first numeric variable if 'DAT' not found
            fn = fieldnames(S);
            got = false;
            for ii = 1:numel(fn)
                if isnumeric(S.(fn{ii}))
                    DAT = S.(fn{ii});
                    got = true;
                    break;
                end
            end
            if ~got, error('smart_binning_1D: No numeric matrix found in MAT file.'); end
        end
    else
        DAT = DAT_in;
    end

    [N, P] = size(DAT);
    K = all_analyte_length;
    if N ~= n_measure * K
        error('smart_binning_1D: size mismatch. N = %d but n_measure*K = %d.', N, n_measure*K);
    end

    % ---- Precompute per-class indices ----
    classIdx = cell(1, K);
    for k = 1:K
        classIdx{k} = ( (k-1)*n_measure + 1 ) : ( k*n_measure );
    end

    % ---- Compute class means/stds per sensor ----
    mu  = zeros(K, P);
    sig = zeros(K, P);
    for j = 1:P
        col = DAT(:, j);
        for k = 1:K
            xk = col(classIdx{k});
            mu(k, j)  = mean(xk, 'omitnan');
            sig(k, j) = std(xk, 0, 'omitnan');
        end
    end

    % ---- Sort classes by mean per sensor; carry stds along ----
    mu_sorted    = zeros(K, P);
    sigma_sorted = zeros(K, P);
    orderMat     = zeros(K, P, 'int32');
    for j = 1:P
        [mu_sorted(:, j), ord] = sort(mu(:, j), 'ascend');
        sigma_sorted(:, j)     = sig(ord, j);
        orderMat(:, j)         = int32(ord(:));
    end

    % ---- Build decision zones per sensor ----
    zoneInfo = repmat(struct('order', [], 'boundaries', [], 'mu', [], 'sigma', []), 1, P);
    Binned   = zeros(N, P, 'int32');

    epsSigma = 1e-12;  % variance floor

    for j = 1:P
        mu_j    = mu_sorted(:, j).';
        sig_j   = max(sigma_sorted(:, j).', sqrt(epsSigma));   % avoid zero std
        ord_j   = double(orderMat(:, j).');

        % Boundaries between adjacent classes in the *sorted* order
        if Decision_boundaries == 1
            % ---- Voronoi (nearest-mean) thresholds: midpoints ----
            b = 0.5 * (mu_j(1:end-1) + mu_j(2:end));
        else
            % ---- QDA thresholds between adjacent 1D Gaussians ----
            b = zeros(1, K-1);
            for t = 1:(K-1)
                b(t) = qda_boundary_1d(mu_j(t), sig_j(t), mu_j(t+1), sig_j(t+1));
            end
        end

        % Store zone info
        zoneInfo(j).order      = ord_j;        % mapping from sorted position -> original class label
        zoneInfo(j).boundaries = b;            % ascending thresholds
        zoneInfo(j).mu         = mu_j;         % sorted means
        zoneInfo(j).sigma      = sig_j;        % sorted stds

        % ---- Assign each sample x to a zone --> original class label ----
        x = DAT(:, j);
        % zone index z in {1..K}: count how many boundaries x exceeds + 1
        % (handles NaNs by default -> gives zone 1; adjust if needed)
        z = ones(N, 1, 'int32');
        if ~isempty(b)
            % For each boundary, increment when x >= boundary
            % Using bsxfun for speed (R2016b+ implicit expansion also fine)
            z = int32(sum(x >= b, 2) + 1);
        end
        % Map sorted zone index -> original class label
        % ord_j(z) gives the final label in {1..K} for this sensor
        Binned(:, j) = int32(ord_j(z).');
    end
end

% ---------- helper: QDA boundary between two 1D Gaussians (equal priors) ----------
function xb = qda_boundary_1d(mu1, s1, mu2, s2)
    % Returns a single boundary between two adjacent classes.
    % If two solutions exist, pick the one between the means; otherwise fallback to midpoint.
    v1 = s1.^2; v2 = s2.^2;
    if abs(v1 - v2) < 1e-15
        xb = 0.5*(mu1 + mu2);             % equal variances -> midpoint
        return;
    end
    % Solve: (x-mu1)^2/(2v1) + 0.5*log v1  =  (x-mu2)^2/(2v2) + 0.5*log v2
    % -> a x^2 + b x + c = 0
    a = (1/(2*v2)) - (1/(2*v1));
    b = (mu1/v1) - (mu2/v2);
    c = (mu2^2)/(2*v2) - (mu1^2)/(2*v1) - 0.5*log(v2/v1);

    if abs(a) < 1e-18
        % Degenerates to linear
        xb = -c / b;
    else
        disc = b*b - 4*a*c;
        if disc < 0
            xb = 0.5*(mu1 + mu2);         % no real root (rare numerically) -> midpoint
        else
            r1 = (-b + sqrt(disc)) / (2*a);
            r2 = (-b - sqrt(disc)) / (2*a);
            % choose the root between the means if possible
            lo = min(mu1, mu2); hi = max(mu1, mu2);
            cand = [r1, r2];
            in  = cand(cand > lo & cand < hi);
            if ~isempty(in)
                % if both are in (pathological), pick the closer to midpoint
                [~, ix] = min(abs(in - 0.5*(mu1+mu2)));
                xb = in(ix);
            else
                % fallback: pick the one closer to the midpoint
                [~, ix] = min(abs(cand - 0.5*(mu1+mu2)));
                xb = cand(ix);
            end
        end
    end
end
