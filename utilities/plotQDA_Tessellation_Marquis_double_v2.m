function plotQDA_Tessellation_Marquis_double_v2(analyte_name_vec, marquis, ...
    classStats_first, classStats, name_title, colorit, x1min, x2min, x1max, x2max)

% 4x4 cm figure; two stacked rectangles; roomy top space for 1-row legend.
figH = figure('Color','w');
set(figH,'Units','centimeters','Position',[2 2 5.9 1.5*5.9/2]);

% ---- layout (normalized) ----
% generous left margin so 'PC2' at 6pt fits;
% shorter panels so the legend has headroom + spare.
left   = 0.2;              % wider left margin for y-label
right  = 0.08;
axW    = 1 - left - right;  % slightly narrower panels
axH    = 0.2;              % shorter panels
gapY   = 0.12;              % gap between panels
botY   = 0.23;              % bottom margin (x-label)
topY   = botY + axH + gapY; % y of top panel -> top edge ~0.78 (room for legend)

axTop = axes('Parent',figH,'Units','normalized','Position',[left topY axW axH]);
axBot = axes('Parent',figH,'Units','normalized','Position',[left botY axW axH]);

% ---- draw (rectangular) ----
[legH_T, legLabels] = drawQDApanel_rect(axTop, classStats_first, analyte_name_vec, ...
                                        marquis, colorit, x1min, x2min, x1max, x2max, 0);
[~, ~]              = drawQDApanel_rect(axBot, classStats,        analyte_name_vec, ...
                                        marquis, colorit, x1min, x2min, x1max, x2max, 1);

% ---- legend (same logic as before; plenty of headroom) ----
h = legend(axTop, legH_T, legLabels, 'Location','NorthOutside','Orientation','horizontal');
h.ItemTokenSize(1) = 3;
set(h,'Units','normalized');

basePos = [0.20 0.8 0.60 0.09];   % unchanged reference box
baseK   = 6;                        % reference count
K       = numel(classStats);
w = min(max(basePos(3) * (K / baseK), 0.25), 0.95);
cx = basePos(1) + basePos(3)/2;
x  = cx - w/2;
h.Position = [x, basePos(2), w, basePos(4)];   % stays high; extra space above panels

% ---- save ----
print(figH, name_title, '-dpng','-r600');
savefig(figH, name_title);
end

% ---------- helpers ----------
function [legH, legLabels] = drawQDApanel_rect(ax, classStats, analyte_name_vec, ...
                                               marquis, colorit, x1min, x2min, x1max, x2max, ylabelon)
axes(ax); hold(ax,'on');
set(ax,'FontSize',6,'LineWidth',0.5,'Box','on', ...
       'LabelFontSizeMultiplier',1);      % keep labels exactly 6 pt

K = numel(classStats); d = 2;
mu = zeros(K,d);  Sigma = zeros(d,d,K);
for k = 1:K
    mu(k,:) = classStats(k).meanXY(:).';
    S       = (classStats(k).covXY + classStats(k).covXY.')/2;
    Sigma(:,:,k) = S;
end
pi_k = ones(K,1)/K;

res = 300;
[x1Grid,x2Grid] = meshgrid(linspace(x1min,x1max,res), linspace(x2min,x2max,res));
[labelGrid, ~]  = qdaPredictFromModels([x1Grid(:), x2Grid(:)], mu, Sigma, pi_k);
labelGrid = reshape(labelGrid, size(x1Grid));

hMap = imagesc(ax, linspace(x1min,x1max,res), linspace(x2min,x2max,res), labelGrid);
set(hMap,'AlphaData',0.30);
axis(ax,'xy'); axis(ax,'tight');              % rectangular (no square/equal)
colormap(ax, colorit(1:K,:));
contour(ax, x1Grid, x2Grid, labelGrid, 'k', 'LineWidth', 0.2);

xlim(ax,[x1min x1max]); ylim(ax,[x2min x2max]);
ylabel(ax,'PC2');
if ylabelon==1, xlabel(ax,'PC1','FontSize',6); end

% Markers
markers = buildMarkers(marquis, K);
legH = gobjects(K,1); legLabels = strings(K,1);
for k=1:K
    h  = scatter(ax, mu(k,1), mu(k,2), 30, 'Marker', markers{k}, ...
        'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',0.2, ...
        'DisplayName', analyzeName(analyte_name_vec,k));
    legH(k) = h; legLabels(k) = string(analyzeName(analyte_name_vec,k));
end
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

function [labels, logPost] = qdaPredictFromModels(Z, mu, Sigma, pi_k)
K = size(mu,1); N = size(Z,1);
logPost = zeros(N,K);
R = cell(K,1); logDet = zeros(1,K);
for k=1:K
    R{k} = chol(Sigma(:,:,k));
    logDet(k) = 2*sum(log(diag(R{k})));
end
for k=1:K
    Y = (Z - mu(k,:)) / R{k};
    logPost(:,k) = -0.5*(sum(Y.^2,2) + logDet(k)) + log(pi_k(k));
end
[~,labels] = max(logPost,[],2);
end
