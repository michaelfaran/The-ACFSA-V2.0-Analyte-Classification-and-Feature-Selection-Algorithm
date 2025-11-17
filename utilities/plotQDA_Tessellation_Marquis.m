function plotQDA_Tessellation_Marquis(analyte_name_vec, marquis, classStats, name_title,colorit,x1min, x2min, x1max, x2max)
% plotQDA_Tessellation_Marquis
%   Build a QDA decision map from Gaussian class models (equal priors),
%   plot the tessellation and ONLY the tile-center markers using the
%   provided marquis string, add legend labels from analyte_name_vec,
%   and save the figure as PNG + FIG using name_title.
%
% Inputs:
%   analyte_name_vec  : cell array of class names for legend, 1xK
%   marquis           : char/string with at least K marker characters
%   classStats        : Kx1 struct with fields:
%                         .meanXY (1x2), .covXY (2x2)
%   name_title        : string/char, used for title and file names
%
% Example:
%   plotQDA_Tessellation_Marquis({'A','B','C','D','E'},'^d<v>d',classStats,'MyPlot')

% ---- extract models
K = numel(classStats);
d = 2;
mu = zeros(K,d);  Sigma = zeros(d,d,K);
for k = 1:K
    mu(k,:)     = classStats(k).meanXY(:).';
    S           = (classStats(k).covXY + classStats(k).covXY.')/2; % symmetrize
    % tiny jitter to ensure PD
    %t           = trace(S)/max(d,1);
    Sigma(:,:,k)= S;
    %+ (1e-9*max(t,1e-12))*eye(d);
end
pi_k = ones(K,1)/K; % equal priors

% ---- global grid bounds (pad around all means using largest std per class)
evMax = zeros(K,1);
for k=1:K, evMax(k) = sqrt(max(eig(Sigma(:,:,k)))); end
% pad = 1;                         % ~5σ truncation
% x1min = min(mu(:,1) - pad*evMax); x1max = max(mu(:,1) + pad*evMax);
% x2min = min(mu(:,2) - pad*evMax); x2max = max(mu(:,2) + pad*evMax);
% bottomX=floor(x1min);
% bottomY=floor(x2min);
% TopX=ceil(x1max);
% TopY=ceil(x2max);
% x2max=TopY;
% x1max=TopX;
% x1min=bottomX;
% x2min=bottomY;

% ---- grid and QDA prediction
res = 300;
[x1Grid,x2Grid] = meshgrid(linspace(x1min,x1max,res), linspace(x2min,x2max,res));
gridPoints = [x1Grid(:), x2Grid(:)];
[labelGrid, logPost] = qdaPredictFromModels(gridPoints, mu, Sigma, pi_k); %#ok<ASGLU>
labelGrid = reshape(labelGrid, size(x1Grid));

% ---- plot tessellation and boundaries
% figure('Color','w');
figgg=figure('Color','w');
bbb=get(figgg,'Position');
new_width=4;
set(figgg, 'Units', 'centimeters', 'Position',[2 2 new_width 1*new_width]);
hold on;
hMap = imagesc(linspace(x1min,x1max,res), linspace(x2min,x2max,res), labelGrid);
set(hMap,'AlphaData',0.30); axis xy; axis equal; axis tight;
% colormap([1 0.8 0.8; 0.8 1 0.8; 0.8 0.8 1; 1 1 0.7; 0.7 1 1]);
colormap(colorit(1:1:K,:));
contour(x1Grid, x2Grid, labelGrid, 'k', 'LineWidth', 0.5);
            ax2=gca;
            grid on;
            axis tight;
xlim([x1min x1max] );
ylim([x2min x2max] );
xlabel('PC1');
ylabel('PC2');
set(gca,'FontSize',6);

% ---- tile centers only (unfilled markers from marquis)
markers = cell(1,K);
fallback = {'o','s','^','v','d','>','<','x','+'};
for k=1:K
    if k <= strlength(string(marquis))
        m = extractBetween(string(marquis), k, k); markers{k} = m{1};
    else
        markers{k} = fallback{min(k,numel(fallback))};
    end
end
colors = lines(K);
legH = gobjects(K,1);

for k=1:K
    mask = (labelGrid == k);
    if any(mask(:))
        % xc = mean(x1Grid(mask)); yc = mean(x2Grid(mask));
        xc =mu(k,1);yc =mu(k,2);
        h  = scatter(xc, yc, 6, 'Marker', markers{k}, ...
             'MarkerEdgeColor', 'k', 'MarkerFaceColor','none', ...
             'LineWidth', 0.5, 'DisplayName', analyzeName(analyte_name_vec,k));
        legH(k) = h;
    else
        % if a tile disappears, still show legend entry
        legH(k) = scatter(nan,nan,90,'Marker',markers{k}, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor','none', ...
            'LineWidth',0.5, 'DisplayName', analyzeName(analyte_name_vec,k));
    end
end

h=legend(legH, 'Location','NorthOutside','Orientation','horizontal');

       ax = gca;
       ax.Position=ax.Position- [0 0 0 0.1];
       % h.Position=[ 0.0435    0.9163    0.9143    0.0488];
       h.ItemTokenSize(1)=3;
       % h.Position= [0.2067    0.8949-0.02   0.6233    0.0861];



% ---- legend (compact, top) ----
% h = legend(ax, legH, legLabels, 'Location','NorthOutside','Orientation','horizontal');
% h.ItemTokenSize(1) = 6;
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
ax2.Position(4)=0.6;
ax2.Box = 'on';         % ensure the box is drawn
ax2.LineWidth = 1; 
% ---- save
grid off
axis square
print(name_title,'-dpng','-r1200');
savefig(name_title);

end % main function

% ---------- helpers ----------
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
% Compute QDA log-posteriors and winning labels given (mu, Sigma, priors)
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

