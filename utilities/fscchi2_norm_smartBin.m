function [ranked,pp] = fscchi2_norm_smartBin(X,Y,varargin)
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

    % internal.stats.checkSupportedNumeric('X',X);
    
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

    % % Initialize binned matrix
    % Xbinned = zeros(size(X),'int32');
    % 
    % % Bin continuous features with per-column normalization
    % for j = 1:D
    %     if iscat(j)
    %         % categorical column
    %         Xbinned(:,j) = classreg.learning.fsutils.indexCategoricals(X(:,j));
    %     else
    %         % continuous column: normalize per column
    %         col = X(:,j);
    %         colMin = min(col);
    %         colMax = max(col);
    % 
    %         if colMax == colMin
    %             % constant column: assign all to bin 1
    %             Xbinned(:,j) = ones(size(col),'int32');
    %         else
    %             % normalize to [0,1] then scale to bins 1..nBins
    %             colNorm = (col - colMin) / (colMax - colMin);
    %             Xbinned(:,j) = min(nBins, max(1, floor(colNorm*nBins)+1));
    %         end
    %     end
    % end
Xbinned=X;
    % Compute p-values using chi-squared test
    p = classreg.learning.fsutils.chi2test(Xbinned,Y,weights,useMissing);
    p = p(:)';

    % Convert to importance scores and rank
    p = -log(p);
    [pp,ranked] = sort(p,'descend'); 