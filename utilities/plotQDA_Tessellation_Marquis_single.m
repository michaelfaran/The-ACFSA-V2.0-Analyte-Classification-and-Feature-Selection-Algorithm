function plotQDA_Tessellation_Marquis_single(analyte_name_vec, marquis, classStats, name_title, colorit, x1min, x2min, x1max, x2max)
% One-square QDA map (quadratic discriminant) with original sizing & legend.

% ---- figure (keep original size) ----
figH = figure('Color','w');
new_width = 3.8; % cm
set(figH,'Units','centimeters','Position',[2 2 new_width 1*new_width]);

% ---- axis layout (slightly smaller square so labels/legend fit) ----
ax = axes('Parent',figH,'Units','normalized','Position',[0.18 0.20 0.64 0.62]);

% ---- draw panel (QDA decision regions) ----
[legH, legLabels] = drawQDAPanel(ax, classStats, analyte_name_vec, marquis, colorit, x1min, x2min, x1max, x2max);

% ---- legend (compact, top; auto-scale width with K) ----
h = legend(ax, legH, legLabels, 'Location','NorthOutside','Orientation','horizontal');
h.ItemTokenSize(1) = 3;
set(h,'Units','normalized');

basePos = [0.20 0.87 0.60 0.09];  % [x y w h]
baseK   = 6;
K       = numel(classStats);
w = basePos(3) * (K / baseK);
w = min(max(w, 0.25), 0.95);
cx = basePos(1) + basePos(3)/2;
x  = cx - w/2;
h.Position = [x, basePos(2), w, basePos(4)];

grid off
% ---- save ----
print(figH, name_title, '-dpng','-r1200');
savefig(figH, name_title);
end

% ---------- helpers ----------
function [legH, legLabels] = drawQDAPanel(ax, classStats, analyte_name_vec, marquis, colorit, x1min, x2min, x1max, x2max)
axes(ax); hold(ax,'on');

% --- Collect class parameters ---
K = numel(classStats);
mu = zeros(K,2);
Sigma = zeros(2,2,K);
pi_k = zeros(K,1);

for k = 1:K
    % mean
    mu_k = classStats(k).meanXY(:).';
    if numel(mu_k) ~= 2
        error('classStats(%d).meanXY must be length 2.', k);
    end
    mu(k,:) = mu_k;

    % covariance (accept several common field names)
    if isfield(classStats(k),'covXY'),    S = classStats(k).covXY;
    elseif isfield(classStats(k),'Sigma'), S = classStats(k).Sigma;
    elseif isfield(classStats(k),'cov'),   S = classStats(k).cov;
    elseif isfield(classStats(k),'covMat'),S = classStats(k).covMat;
    else, error('Provide 2x2 covariance for class %d (fields: covXY/Sigma/cov/covMat).', k);
    end
    if ~isequal(size(S),[2 2])
        error('Covariance for class %d must be 2x2.', k);
    end
    Sigma(:,:,k) = S;

    % prior (optional; default equal)
    if isfield(classStats(k),'prior') && ~isempty(classStats(k).prior)
        pi_k(k) = classStats(k).prior;
    else
        pi_k(k) = 1; % will normalize later
    end
end
pi_k = pi_k / sum(pi_k);

% --- Grid for classification ---
res = 300;
[x1Grid, x2Grid] = meshgrid(linspace(x1min,x1max,res), linspace(x2min,x2max,res));
GX = x1Grid(:); GY = x2Grid(:);
X = [GX GY];                        % N x 2
N = size(X,1);

% --- QDA discriminants: g_k(x) = -0.5 log|Σ_k| - 0.5 (x-μ_k)^T Σ_k^{-1} (x-μ_k) + log π_k ---
G = -inf(N,K);
for k = 1:K
    S = Sigma(:,:,k);

    % Regularize if needed for numerical stability
    S = regularize2x2(S);

    % Chol for logdet & solving
    [R, p] = chol(S);
    if p ~= 0
        % In rare case still not PD, add small ridge and retry
        S = S + 1e-8*eye(2);
        [R, p] = chol(S);
        if p ~= 0
            % last resort: symmetrize
            S = (S+S.')/2;
            [R,p] = chol(S + 1e-8*eye(2));
            if p~=0, error('Covariance for class %d is not PD even after regularization.', k); end
        end
    end

    logdetS = 2*sum(log(diag(R)));      % log determinant via chol
    Xm = bsxfun(@minus, X, mu(k,:));    % N x 2
    % Mahalanobis term: (x-μ)ᵗ Σ^{-1} (x-μ) = sum( (Xm / S) .* Xm, 2 )
    Mah = sum((Xm / S) .* Xm, 2);

    G(:,k) = -0.5*(Mah + logdetS) + log(max(pi_k(k), realmin));
end

[~, labels] = max(G, [], 2);
labelGrid = reshape(labels, size(x1Grid));

% --- Plot decision regions & boundaries ---
hMap = imagesc(ax, linspace(x1min,x1max,res), linspace(x2min,x2max,res), labelGrid);
set(hMap,'AlphaData',0.30); axis(ax,'xy'); axis(ax,'equal'); axis(ax,'tight');
if size(colorit,1) >= K
    colormap(ax, colorit(1:K,:));
end
contour(ax, x1Grid, x2Grid, labelGrid, 'k', 'LineWidth', 0.2);

grid(ax,'on');
xlim(ax,[x1min x1max]); ylim(ax,[x2min x2max]);
xlabel(ax,'PC1'); ylabel(ax,'PC2');
set(ax,'FontSize',6,'LineWidth',0.2,'Box','on');
axis(ax,'square');

% --- Class mean markers (same as Voronoi version) ---
markers = buildMarkers(marquis, K);
legH = gobjects(K,1);
legLabels = strings(K,1);
for k=1:K
    h  = scatter(ax, mu(k,1), mu(k,2), 16, 'Marker', markers{k}, ...
        'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',0.1, ...
        'DisplayName', analyzeName(analyte_name_vec,k));
    legH(k) = h;
    legLabels(k) = string(analyzeName(analyte_name_vec,k));
end
end

function Sreg = regularize2x2(S)
% Gentle, data-scaled ridge to avoid singularities with tiny sample sizes
S = (S+S')/2;                   % enforce symmetry
t = trace(S);
if ~isfinite(t) || t<=0, t = 1; end
lambda = 1e-8 * t/2 + 1e-12;
Sreg = S + lambda*eye(2);
end

function markers = buildMarkers(marquis, K)
markers = cell(1,K);
fallback = {'o','s','^','v','d','>','<','x','+'};
for k=1:K
    if k <= strlength(string(marquis))
        m = extractBetween(string(marquis), k, k); markers{k} = m{1};
    else
        markers{k} = fallback{min(k,numel(fallback))};
    end
end
end

function name = analyzeName(analyte_name_vec, k)
    if iscell(analyte_name_vec) && numel(analyte_name_vec) >= k
        name = analyte_name_vec{k};
    elseif isstring(analyte_name_vec) && numel(analyte_name_vec) >= k
        name = analyte_name_vec(k);
    else
        name = sprintf('Class %d',k);
    end
end
