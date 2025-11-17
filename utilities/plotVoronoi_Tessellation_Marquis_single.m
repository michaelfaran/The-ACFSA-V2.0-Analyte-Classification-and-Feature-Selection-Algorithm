function plotVoronoi_Tessellation_Marquis_single(analyte_name_vec, marquis, classStats, name_title, colorit, x1min, x2min, x1max, x2max)
% One square Voronoi map (nearest-mean) with original sizing & legend.

% ---- figure (keep original size) ----
figH = figure('Color','w');
new_width = 3.8; % cm
set(figH,'Units','centimeters','Position',[2 2 new_width 1*new_width]);

% ---- axis layout (slightly smaller square so labels/legend fit) ----
ax = axes('Parent',figH,'Units','normalized','Position',[0.18 0.20 0.64 0.62]);

% ---- draw panel ----
[legH, legLabels] = drawVoronoiPanel(ax, classStats, analyte_name_vec, marquis, colorit, x1min, x2min, x1max, x2max);

% ---- legend (compact, top) ----
h = legend(ax, legH, legLabels, 'Location','NorthOutside','Orientation','horizontal');
h.ItemTokenSize(1) = 3;
set(h,'Units','normalized');
% h.Position = [0.20 0.87 0.60 0.09];

% Base legend box (your original numbers)
basePos = [0.20 0.87 0.60 0.09];   % [x y w h] in normalized units
baseK   = 6;                                % reference count
K       = numel(classStats);

% Scale width proportional to item count, clamp to sane bounds
w = basePos(3) * (K / baseK);
w = min(max(w, 0.25), 0.95);               % clamp between 0.25 and 0.95

% Keep legend horizontally centered around original center
cx = basePos(1) + basePos(3)/2;
x  = cx - w/2;

% Apply new position (same y & height)
h.Units    = 'normalized';
h.Position = [x, basePos(2), w, basePos(4)];


grid off
% ---- save ----
print(figH, name_title, '-dpng','-r1200');
savefig(figH, name_title);
end

% ---------- helpers ----------
function [legH, legLabels] = drawVoronoiPanel(ax, classStats, analyte_name_vec, marquis, colorit, x1min, x2min, x1max, x2max)
axes(ax); hold(ax,'on');

% Means (Voronoi sites)
K = numel(classStats);
mu = zeros(K,2);
for k = 1:K
    mu(k,:) = classStats(k).meanXY(:).';
end

% Grid & nearest-mean labeling
res = 300;
[x1Grid, x2Grid] = meshgrid(linspace(x1min,x1max,res), linspace(x2min,x2max,res));
GX = x1Grid(:); GY = x2Grid(:);

D2 = zeros(numel(GX), K);
for k = 1:K
    dx = GX - mu(k,1);
    dy = GY - mu(k,2);
    D2(:,k) = dx.*dx + dy.*dy;
end
[~, labels] = min(D2, [], 2);
labelGrid = reshape(labels, size(x1Grid));

% Plot tessellation
hMap = imagesc(ax, linspace(x1min,x1max,res), linspace(x2min,x2max,res), labelGrid);
set(hMap,'AlphaData',0.30); axis(ax,'xy'); axis(ax,'equal'); axis(ax,'tight');
colormap(ax, colorit(1:K,:));
contour(ax, x1Grid, x2Grid, labelGrid, 'k', 'LineWidth', 0.2);

grid(ax,'on');
xlim(ax,[x1min x1max]); ylim(ax,[x2min x2max]);
xlabel(ax,'PC1'); ylabel(ax,'PC2');
set(ax,'FontSize',6,'LineWidth',0.2,'Box','on');
axis(ax,'square');

% Markers at class means
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
