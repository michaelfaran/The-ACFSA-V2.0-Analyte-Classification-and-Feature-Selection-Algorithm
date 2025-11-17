function plotQDA_Tessellation_Marquis_double(analyte_name_vec, marquis, classStats_first, classStats, name_title, colorit, x1min, x2min, x1max, x2max)
% Two square QDA maps side-by-side with adjusted sizing so labels and legend fit.

% ---- figure (keep original size) ----
figH = figure('Color','w');
new_width = 4; % cm
set(figH,'Units','centimeters','Position',[2 2 new_width 1*new_width]);

% ---- axes layout (two squares, slightly smaller to leave room) ----
axL = axes('Parent',figH,'Units','normalized','Position',[0.12 0.20 0.33 0.62]); % left
axR = axes('Parent',figH,'Units','normalized','Position',[0.55 0.20 0.33 0.62]); % right
% ---- draw both panels ----
[legH_L, legLabels] = drawQDApanel(axL, classStats_first, analyte_name_vec, marquis, colorit, x1min, x2min, x1max, x2max,1);
grid off;
[legH_R, ~]        = drawQDApanel(axR, classStats,        analyte_name_vec, marquis, colorit, x1min, x2min, x1max, x2max,1);
grid off;
% % ---- legend (same style as your original, but now fits above axes) ----
 h = legend(axR, legH_R, legLabels, 'Location','NorthOutside','Orientation','horizontal');
% h.ItemTokenSize(1) = 6;
% set(h,'Units','normalized');
% h.Position = [0.20 0.87 0.60 0.09]; % shifted a bit up to avoid overlap


% ---- legend (compact, top) ----
% h = legend(ax, legH, legLabels, 'Location','NorthOutside','Orientation','horizontal');
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

% ---- save ----
print(figH, name_title, '-dpng','-r600');
savefig(figH, name_title);

end % main

% ---------- helpers ----------
function [legH, legLabels] = drawQDApanel(ax, classStats, analyte_name_vec, marquis, colorit, x1min, x2min, x1max, x2max,ylabelon)
axes(ax); hold(ax,'on');

% Extract models
K = numel(classStats);
d = 2;
mu = zeros(K,d);  Sigma = zeros(d,d,K);
for k = 1:K
    mu(k,:)      = classStats(k).meanXY(:).';
    S            = (classStats(k).covXY + classStats(k).covXY.')/2; % symmetrize
    Sigma(:,:,k) = S;
end
pi_k = ones(K,1)/K;

% Grid + QDA prediction
res = 300;
[x1Grid,x2Grid] = meshgrid(linspace(x1min,x1max,res), linspace(x2min,x2max,res));
gridPoints = [x1Grid(:), x2Grid(:)];
[labelGrid, ~] = qdaPredictFromModels(gridPoints, mu, Sigma, pi_k);
labelGrid = reshape(labelGrid, size(x1Grid));

% Tessellation
hMap = imagesc(ax, linspace(x1min,x1max,res), linspace(x2min,x2max,res), labelGrid);
set(hMap,'AlphaData',0.30); axis(ax,'xy'); axis(ax,'equal'); axis(ax,'tight');
colormap(ax, colorit(1:K,:));
contour(ax, x1Grid, x2Grid, labelGrid, 'k', 'LineWidth', 0.5);

grid(ax,'on');
xlim(ax,[x1min x1max]); ylim(ax,[x2min x2max]);
xlabel(ax,'PC1'); 
if ylabelon==1
ylabel(ax,'PC2');
end
set(ax,'FontSize',6,'LineWidth',1,'Box','on');
axis(ax,'square');

% Markers
markers = buildMarkers(marquis, K);
legH = gobjects(K,1);
legLabels = strings(K,1);
for k=1:K
    xc = mu(k,1); yc = mu(k,2);
    h  = scatter(ax, xc, yc, 16, 'Marker', markers{k}, ...
        'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',0.5, ...
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

function [labels, logPost] = qdaPredictFromModels(Z, mu, Sigma, pi_k)
K = size(mu,1); N = size(Z,1);
logPost = zeros(N,K);
R = cell(K,1); logDet = zeros(1,K);
for k=1:K
    R{k} = chol(Sigma(:,:,k));
    logDet(k) = 2*sum(log(diag(R{k})));
end
for k=1:K
    XC = Z - mu(k,:);
    Y  = XC / R{k};
    Q  = sum(Y.^2,2);
    logPost(:,k) = -0.5*(Q + logDet(k)) + log(pi_k(k));
end
[~,labels] = max(logPost,[],2);
end
