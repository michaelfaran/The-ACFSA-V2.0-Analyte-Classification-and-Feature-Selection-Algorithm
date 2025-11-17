function [muMat, varMat, new_DAT] = analyte_stats_resample(DAT, n_measure, n_resample, rng_seed, STD_buff_new)
% ANALYTE_STATS_RESAMPLE
% DAT:  (n_analytes*n_measure) x n_sensors; rows grouped by analyte (blocks of n_measure)
% Returns:
%   muMat   : n_analytes x n_sensors  (per-analyte, per-sensor means)
%   varMat  : n_analytes x n_sensors  (per-analyte, per-sensor variances; unbiased)
%   new_DAT : (n_analytes*n_resample) x n_sensors synthetic data
%
% n_resample    : number of synthetic samples per analyte (default 20)
% rng_seed      : optional RNG seed
% STD_buff_new  : factor to scale STD when sampling new_DAT (default 1)
%                 - scalar: same factor for all (analyte, sensor)
%                 - 1×n_sensors vector: per-sensor factor (broadcast over analytes)
%                 - n_analytes×n_sensors matrix: per (analyte, sensor) factor

    if nargin < 3 || isempty(n_resample), n_resample = 20; end
    if nargin >= 4 && ~isempty(rng_seed), rng(rng_seed); end
    if nargin < 5 || isempty(STD_buff_new), STD_buff_new = 1; end

    [n_rows, n_sensors] = size(DAT);
    assert(mod(n_rows, n_measure) == 0, 'Rows must be a multiple of n_measure.');
    n_analytes = n_rows / n_measure;

    % reshape to (n_measure x n_sensors x n_analytes)
    data3 = permute(reshape(DAT.', n_sensors, n_measure, n_analytes), [2 1 3]);

    % per-analyte, per-sensor stats
    muMat  = squeeze(mean(data3, 1, 'omitnan')).';    % n_analytes x n_sensors
    varMat = squeeze(var( data3, 0, 1, 'omitnan')).'; % unbiased (N-1)

    % build STD buffer matrix
    if isscalar(STD_buff_new)
        stdBuff = repmat(STD_buff_new, n_analytes, n_sensors);
    elseif isvector(STD_buff_new) && numel(STD_buff_new) == n_sensors
        stdBuff = repmat(reshape(STD_buff_new,1,[]), n_analytes, 1);
    elseif isequal(size(STD_buff_new), [n_analytes, n_sensors])
        stdBuff = STD_buff_new;
    else
        error('STD_buff_new must be scalar, 1xN_sensors, or N_analytesxN_sensors.');
    end
    if any(stdBuff(:) < 0)
        error('STD_buff_new must be nonnegative.');
    end

    % synthesize new data with inflated/shrunk STD
    new_DAT = nan(n_analytes * n_resample, n_sensors);
    for a = 1:n_analytes
        rows = (a-1)*n_resample + (1:n_resample);
        for s = 1:n_sensors
            mu = muMat(a,s);
            v  = varMat(a,s);
            sf = stdBuff(a,s);
            if isnan(mu) || isnan(v) || v < eps || sf==0
                new_DAT(rows, s) = mu; % constant if variance ~0 or scaling is zero
            else
                new_DAT(rows, s) = mu + (sqrt(v)*sf) * randn(n_resample,1);
            end
        end
    end
end
