function [ranked,p] = fscchi2_norm_UP(X, Y, classStats, varargin)
%FSCCHI2_NORM Univariate feature selection for classification using Chi-2 test
%   Modified version with per-column min-max normalization for continuous features.

    if nargin > 1
        Y = convertStringsToChars(Y);
    end

    if nargin > 2
        [varargin{:}] = convertStringsToChars(varargin{:});
    end

    args = {'NumBins','UseMissing'};
    defaults = {10,false};
    [nBins,useMissing,~,otherArgs] = internal.stats.parseArgs(args,defaults,varargin{:});
    
    % Error conditions
    if ~isnumeric(nBins) || ~isscalar(nBins) ...
            || ~isreal(nBins) || nBins~=round(nBins) ...
            || nBins<=0 || isnan(nBins) || isinf(nBins)
        error(message('stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadNumBins'));
    end

    useMissing = internal.stats.parseOnOff(useMissing,'UseMissing');
    
    [X,Y,weights,dataSummary] = ...
      classreg.learning.classif.FullClassificationModel.prepareData(...
      X,Y,otherArgs{:},'OrdinalIsCategorical',false);
  
    Y = int32(grp2idx(Y));
    D = size(X,2);

    internal.stats.checkSupportedNumeric('X',X);
    
    % Edge case: one class only or only one bin
    if max(Y)==1 || nBins == 1
        ranked = 1:D;
        p = zeros(1,D);
        return
    end

    % Get indices of categorical features
    catpreds = dataSummary.CategoricalPredictors;
    iscat = false(1,D);
    iscat(catpreds) = true;

    % Initialize binned matrix
    Xbinned = zeros(size(X),'int32');

    % Bin continuous features with per-column normalization
    for j = 1:D
        if iscat(j)
            % categorical column
            Xbinned(:,j) = classreg.learning.fsutils.indexCategoricals(X(:,j));
        else
            % continuous column: normalize per column
            col = X(:,j);
            colMin = min(col);
            colMax = max(col);
            
            if colMax == colMin
                % constant column: assign all to bin 1
                Xbinned(:,j) = ones(size(col),'int32');
            else
                % normalize to [0,1] then scale to bins 1..nBins
                colNorm = (col - colMin) / (colMax - colMin);
                Xbinned(:,j) = min(nBins, max(1, floor(colNorm*nBins)+1));
            end
        end
    end

    % Compute p-values using chi-squared test
    p = classreg.learning.fsutils.chi2test(Xbinned,Y,weights,useMissing);
    p = p(:)';

    % Convert to importance scores and rank
    p = -log(p);
    [pp2,ranked2] = sort(p,'descend'); 


    % -------------------------
    % Pairwise projected Mahalanobis distance
    % -------------------------
    classes = unique(Y);
    nClasses = numel(classes);
    Wa=0.5*(nClasses^2-nClasses);
    final_pairwise_scores = zeros(1,length(p)); % final feature score
    
    for c1 = 1:nClasses-1
        for c2 = c1+1:nClasses

             % Extract samples for this pair
            idxPair = Y == c1 | Y == c2;
            Xpair = Xbinned(idxPair,:);
            Ypair = Y(idxPair);
            
            % Map class labels to 1 and 2
            Ypair_relabeled = zeros(size(Ypair));
            Ypair_relabeled(Ypair == c1) = 1;
            Ypair_relabeled(Ypair == c2) = 2;
            weights=ones(size(Ypair)).*(1./length(Ypair));
            % Compute chi² per feature for this pair
            Ypair_relabeled = int32(Ypair_relabeled);  % or categorical(Ypair_relabeled)
            chi_score_pair = classreg.learning.fsutils.chi2test(Xpair, Ypair_relabeled,weights,useMissing);
            chi_score_pair = chi_score_pair(:)'; % 1xD
             %chi_score_pair = -log(chi_score_pair);
            unique_vec_count=zeros(1,size(chi_score_pair,2));
            for kc=1:1:size(chi_score_pair,2)
                U=unique(Xpair(:,kc));
                unique_vec_count(kc)=length(U);
            end
            chi2_val = chi2inv(1 - chi_score_pair, (2-1).*unique_vec_count);%assuming euqal to numebr of DoF of chi sqaured here
            % Vector between class means in PCA space
            mean1 = classStats(c1).meanXY;
            mean2 = classStats(c2).meanXY;
            v = mean2 - mean1;
            v_norm = v / norm(v);

            % Project each class covariance along vector v
            sigma1_sq = v_norm * classStats(c1).covXY * v_norm';
            sigma2_sq = v_norm * classStats(c2).covXY * v_norm';

            % Combined standard deviation
            sigma_total = sqrt(0.5*(sigma1_sq + sigma2_sq));%discrminality index

            % Projected Mahalanobis distance
            maha_dist = norm(v) / sigma_total;

            % Add contribution to final feature scores
            momo2=(1./((maha_dist)+10^-5)).*chi2_val;
            final_pairwise_scores=final_pairwise_scores+momo2;
            % final_pairwise_scores = final_pairwise_scores + 1 ./ (maha_dist.*chi2_val+1*10^-5.*maha_dist);
        end
    end

    % Rank features by final combined score
    % [nemo, ranked] = sort(1./final_pairwise_scores, 'descend');
    [nemo, ranked] = sort(final_pairwise_scores, 'descend');

end
