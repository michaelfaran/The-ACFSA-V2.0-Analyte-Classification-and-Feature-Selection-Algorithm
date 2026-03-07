function [ranked,pScore] = fscchi2_norm_UP_smartBin(X, Y, classStats, varargin)
%FSCCHI2_NORM_UP_SMARTBIN
% Univariate feature scoring using chi^2 (smart-binned inputs) + pairwise
% covariance-projected separation weighting (uses classStats mean/cov).
%
% IMPORTANT:
%   - X is assumed to already be integer "bin labels" per feature (smart binning),
%     i.e., X(:,j) in {1..K} or at least discrete labels.
%   - classStats(k).meanXY and classStats(k).covXY should already include your
%     regularization (LW shrinkage + eigen-floor ridge), and this function will
%     respect it.
%
% Outputs:
%   ranked  : feature indices sorted best->worst (descending score)
%   pScore  : 1xD score vector (default: -log(pvals_all_classes))

    if nargin > 1
        Y = convertStringsToChars(Y);
    end
    if nargin > 3
        [varargin{:}] = convertStringsToChars(varargin{:});
    end

    args = {'NumBins','UseMissing'};
    defaults = {10,false};
    [nBins,useMissing,~,otherArgs] = internal.stats.parseArgs(args,defaults,varargin{:});

    if ~isnumeric(nBins) || ~isscalar(nBins) || ~isreal(nBins) || nBins~=round(nBins) || nBins<=0 || isnan(nBins) || isinf(nBins)
        error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadNumBins'));
    end
    useMissing = internal.stats.parseOnOff(useMissing,'UseMissing');

    % PrepareData returns weights; DO NOT overwrite them later
    [X,Y,weights,~] = classreg.learning.classif.FullClassificationModel.prepareData( ...
        X, Y, otherArgs{:}, 'OrdinalIsCategorical', false);

    Y = int32(grp2idx(Y));
    D = size(X,2);

    % Edge cases
    if max(Y)==1 || nBins==1
        ranked = 1:D;
        pScore = zeros(1,D);
        return
    end

    % X is already binned (smart-bin labels), keep name for clarity
    Xbinned = X;

    % ---------- global chi^2 p-values (all classes) ----------
    pvals = classreg.learning.fsutils.chi2test(Xbinned, Y, weights, useMissing);
    pvals = pvals(:)';

    % Clamp p-values to avoid -log(0) and chi2inv(1, dof) issues
    pvals = max(min(pvals, 1 - realmin), realmin);

    % Score = -log(p)
    pScore = -log(pvals);

    % ---------- pairwise weighting using classStats ----------
    classes  = unique(Y);
    nClasses = numel(classes);

    final_pairwise_scores = zeros(1, D);

    epsSmall = 1e-12;

    for c1 = 1:(nClasses-1)
        for c2 = (c1+1):nClasses

            % Extract samples for this pair
            idxPair = (Y == classes(c1)) | (Y == classes(c2));
            Xpair   = Xbinned(idxPair,:);
            Ypair   = Y(idxPair);
            wpair   = weights(idxPair);

            % Relabel to {1,2}
            Ypair2 = zeros(size(Ypair), 'int32');
            Ypair2(Ypair == classes(c1)) = 1;
            Ypair2(Ypair == classes(c2)) = 2;

            % Pairwise chi^2 p-values per feature
            p_pair = classreg.learning.fsutils.chi2test(Xpair, Ypair2, wpair, useMissing);
            p_pair = p_pair(:)';

            % Clamp
            p_pair = max(min(p_pair, 1 - realmin), realmin);

            % Convert p->chi2 statistic using dof ≈ (#uniqueBins-1) for 2-class
            chi2_val = zeros(1,D);
            for kc = 1:D
                u = unique(Xpair(:,kc));
                nU = numel(u);
                dof = max(1, nU - 1); % since (r-1)(c-1) and c=2
                chi2_val(kc) = chi2inv(1 - p_pair(kc), dof);
            end

            % Separation along mean-difference direction in PCA space
            mean1 = classStats(c1).meanXY;
            mean2 = classStats(c2).meanXY;
            v = mean2 - mean1;

            nv = norm(v);
            if ~(isfinite(nv) && nv > epsSmall)
                % If means coincide, this pair gives no directional separation
                continue;
            end
            v_norm = v / nv;

            % Project cov along v (robust)
            C1 = (classStats(c1).covXY + classStats(c1).covXY')/2;
            C2 = (classStats(c2).covXY + classStats(c2).covXY')/2;

            sigma1_sq = v_norm * C1 * v_norm';
            sigma2_sq = v_norm * C2 * v_norm';

            sigma1_sq = max(sigma1_sq, epsSmall);
            sigma2_sq = max(sigma2_sq, epsSmall);

            sigma_total = sqrt(0.5*(sigma1_sq + sigma2_sq) + epsSmall);

            maha_dist = nv / sigma_total;   % larger = better separation
            maha_dist = max(maha_dist, epsSmall);

            % Your weighting form (kept), with safety
            w_maha = 1 / (sqrt(maha_dist) + 1e-8);

            % Accumulate feature scores
            final_pairwise_scores = final_pairwise_scores + w_maha .* chi2_val;
        end
    end

    % Final ranking: descending combined score
    [~, ranked] = sort(final_pairwise_scores, 'descend');
end
