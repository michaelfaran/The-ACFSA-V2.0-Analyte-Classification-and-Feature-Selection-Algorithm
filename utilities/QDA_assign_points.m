function [labels, dbg] = QDA_assign_points(X2, classStats, varargin)
% QDA_assign_points  Classify 2D points using QDA built from classStats
%
%   labels = QDA_assign_points(X2, classStats)
%   labels = QDA_assign_points(X2, classStats, 'Res',300,'Trunc',5)
%   labels = QDA_assign_points(..., 'Plot',true, 'PlotCfg',cfg)
%
% Inputs:
%   X2         : N x 2 array of query points
%   classStats : Kx1 struct with fields:
%                  .meanXY  (1x2)
%                  .covXY   (2x2)  -- symmetric positive (semi)definite
% Name-Value (optional):
%   'Res'      : grid resolution (kept for signature compatibility; default 300)
%   'Trunc'    : truncation in sigmas (kept for signature compatibility; default 5)
%   'Plot'     : if true, draws QDA tessellation + overlays the query points (default false)
%   'PlotCfg'  : struct with fields required by plotQDA_Tessellation_Marquis_single:
%                .analyte_name_vec
%                .marquis
%                .name_title
%                .colorit
%                .x1min .x2min .x1max .x2max
%                Optional extras:
%                .PointSize (default 18)
%                .PointEdgeColor (default 'k')
%                .PointFaceColor (default 'none')
%                .Annotate (default false)   % write predicted class index near each point
%
% Output:
%   labels     : N x 1 vector of class indices in {1,...,K}
%   dbg        : (optional) debug struct with logPost, mu, Sigma, etc.

% ---- parse options
p = inputParser;
p.addParameter('Res', 300, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Trunc', 5, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Plot', false, @(b)islogical(b)&&isscalar(b));
p.addParameter('PlotCfg', struct(), @(s)isstruct(s)&&isscalar(s));
p.parse(varargin{:});
opt = p.Results;

% ---- checks
if ~isnumeric(X2) || size(X2,2) ~= 2
    error('X2 must be N x 2 numeric.');
end

% ---- extract class parameters
K = numel(classStats); d = 2;
mu = zeros(K,d);  Sigma = zeros(d,d,K);
for k = 1:K
    mu(k,:) = classStats(k).meanXY(:).';
    S = (classStats(k).covXY + classStats(k).covXY.')/2;  % symmetrize

    % tiny jitter to ensure PD for chol:
    t = trace(S)/max(d,1);
    Sigma(:,:,k) = S + (1e-12 + 1e-9*max(t,1e-12))*eye(d);
end
pi_k = ones(K,1)/K;  % equal priors

% ---- precompute Cholesky + log|Sigma|
R = cell(K,1); logDet = zeros(1,K);
for k = 1:K
    try
        R{k} = chol(Sigma(:,:,k),'upper');  % Sigma = R' * R
    catch
        % if still not PD, nudge a bit more
        bump = 1e-6*trace(Sigma(:,:,k))/d;
        R{k} = chol(Sigma(:,:,k)+bump*eye(d),'upper');
        Sigma(:,:,k) = Sigma(:,:,k)+bump*eye(d);
    end
    logDet(k) = 2*sum(log(diag(R{k})));
end

% ---- compute QDA log-posteriors and pick argmax
N = size(X2,1);
logPost = zeros(N,K);
for k = 1:K
    XC = X2 - mu(k,:);
    % Mahalanobis via triangular solve: (X2-mu)/R, then squared norm
    Y = XC / R{k};
    Q = sum(Y.^2,2);
    logPost(:,k) = -0.5*(Q + logDet(k)) + log(pi_k(k));
end
[~, labels] = max(logPost, [], 2);

% ---- optional plotting sanity check
if opt.Plot
    cfg = opt.PlotCfg;

    % minimal validation (so it fails loudly if something is missing)
    required = {'analyte_name_vec','marquis','name_title','colorit','x1min','x2min','x1max','x2max'};
    for i = 1:numel(required)
        if ~isfield(cfg, required{i})
            error('PlotCfg is missing field "%s".', required{i});
        end
    end

    % Draw tessellation (this creates a new figure)
    plotQDA_Tessellation_Marquis_single( ...
        cfg.analyte_name_vec, cfg.marquis, classStats, cfg.name_title, cfg.colorit, ...
        cfg.x1min, cfg.x2min, cfg.x1max, cfg.x2max);

    ax = gca; hold(ax,'on');

    % point style defaults
    if ~isfield(cfg,'PointSize'),      cfg.PointSize = 18; end
    if ~isfield(cfg,'PointEdgeColor'), cfg.PointEdgeColor = 'k'; end
    if ~isfield(cfg,'PointFaceColor'), cfg.PointFaceColor = 'none'; end
    if ~isfield(cfg,'Annotate'),       cfg.Annotate = false; end

    % overlay points
    scatter(ax, X2(:,1), X2(:,2), cfg.PointSize, ...
        'o', 'LineWidth',0.6, ...
        'MarkerEdgeColor', cfg.PointEdgeColor, ...
        'MarkerFaceColor', cfg.PointFaceColor);

    % optionally annotate predictions
    if cfg.Annotate
        for i = 1:N
            text(ax, X2(i,1), X2(i,2), sprintf('  %d', labels(i)), ...
                'FontSize',6, 'Color','k', 'HorizontalAlignment','left', 'VerticalAlignment','middle');
        end
    end

    hold(ax,'off');
end

% ---- optional debug output
if nargout > 1
    dbg = struct();
    dbg.logPost = logPost;
    dbg.mu      = mu;
    dbg.Sigma   = Sigma;
    dbg.pi_k    = pi_k;
    dbg.logDet  = logDet;
else
    dbg = [];
end
end