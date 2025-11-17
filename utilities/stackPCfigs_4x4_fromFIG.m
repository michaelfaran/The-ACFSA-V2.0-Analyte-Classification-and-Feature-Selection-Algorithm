function figOut = stackPCfigs_4x4_fromFIG(figfile_top, figfile_bottom, out_name)
% Create a 4x4 cm combined figure from two existing .fig files
% (top and bottom panels). Keeps a single one-row legend from the TOP file.
%
% Inputs:
%   figfile_top     - path or handle to first figure (.fig)  -> top panel
%   figfile_bottom  - path or handle to second figure (.fig) -> bottom panel
%   out_name        - filename prefix (no extension) for saving
%
% Output:
%   figOut - handle to the combined figure

    %----- Open (invisible) source figures if paths are given -----
    needClose1 = false; needClose2 = false;
    if ischar(figfile_top) || isstring(figfile_top)
        src1 = openfig(figfile_top, 'invisible'); needClose1 = true;
    else
        src1 = figure(figfile_top); set(src1,'Visible','off');
    end
    if ischar(figfile_bottom) || isstring(figfile_bottom)
        src2 = openfig(figfile_bottom, 'invisible'); needClose2 = true;
    else
        src2 = figure(figfile_bottom); set(src2,'Visible','off');
    end

    %----- Find main axes in each source figure -----
    ax1 = getMainAxes(src1);
    ax2 = getMainAxes(src2);

    %----- Build destination canvas (4 x 4 cm) -----
    figOut = figure('Color','w');
    set(figOut,'Units','centimeters','Position',[2 2 5.9 1.5*5.9/2]);

    % % Layout (matches your v2): narrow + short rectangles, room for legend
    % left   = 0.20;          % wider left margin for 'PC2' label at 6pt
    % right  = 0.08;
    % axW    = 1 - left - right;
    % axH    = 0.23;
    % gapY   = 0.10;
    % botY   = 0.20;
    % topY   = botY + axH + gapY;
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


    axTop = axes('Parent',figOut, 'Units','normalized', 'Position',[left topY axW axH]);
    axBot = axes('Parent',figOut, 'Units','normalized', 'Position',[left botY axW axH]);

    %----- Copy content (children) from each source axes -----
    copyobj(allchild(ax1), axTop);
    copyobj(allchild(ax2), axBot);

    % Copy key axes properties (limits, labels, box/line widths/font sizes)
    cloneAxesProps(ax1, axTop);
    cloneAxesProps(ax2, axBot);

    % Force labels font to 6pt as in your layout
    % set(get(axTop,'XLabel'),'FontSize',6);  
    axTop.XLabel.Visible = 'off';    
    set(get(axTop,'YLabel'),'FontSize',6);
    set(get(axBot,'XLabel'),'FontSize',6);  set(get(axBot,'YLabel'),'FontSize',6);
    set(axTop,'FontSize',6,'LineWidth',0.5,'Box','on');
    set(axBot,'FontSize',6,'LineWidth',0.5,'Box','on');

    %----- Legend: keep ONLY ONE (from TOP source fig) -----
    legSrc = findobj(src1,'Type','Legend');
    if ~isempty(legSrc)
        legSrc = legSrc(1);  % take the first legend
        labels = get(legSrc,'String');           % label text cell array
        pcs    = get(legSrc,'PlotChildren');     % handles used by the legend
        % Legend uses reverse order; flip to display order
        %pcs = flipud(pcs);
        % Recreate "ghost" handles in the TOP destination axes to match tokens
        L = numel(pcs);
        ph = gobjects(L,1);
        for i = 1:L
            ph(i) = copyobj(pcs(i), axTop);      % copy legend token style
            try
                % Ensure they don't render on data (NaN keeps them invisible in plot)
                if isprop(ph(i),'XData'); set(ph(i),'XData',NaN,'YData',NaN); end
                if isprop(ph(i),'ZData'); set(ph(i),'ZData',NaN); end
            end
        end
        % Make one-row legend with your previous positioning logic
        h = legend(axTop, ph, labels, 'Location','NorthOutside', 'Orientation','horizontal');
        h.ItemTokenSize(1) = 3;
        h.Units = 'normalized';

        basePos = [0.20 0.80 0.60 0.09];  % your prior reference box
        baseK   = 6;                      % reference count
        K       = numel(labels);
        w = min(max(basePos(3) * (K / baseK), 0.25), 0.95);
        cx = basePos(1) + basePos(3)/2;
        x  = cx - w/2;
        h.Position = [x, basePos(2), w, basePos(4)];
    end

    %----- Save outputs -----
    if nargin >= 3 && ~isempty(out_name)
        print(figOut, out_name, '-dpng','-r600');
        savefig(figOut, out_name);
    end

    %----- Clean up source figures if we opened them -----
    if needClose1, close(src1); end
    if needClose2, close(src2); end
end

% ========= Helpers =========
function ax = getMainAxes(figH)
    % Return the "main" axes in a figure (ignore overlay/tagged axes)
    axs = findobj(figH,'Type','axes');
    if isempty(axs), error('No axes found in source figure.'); end
    % Drop overlay axes if present
    axs(strcmp(get(axs,'Tag'),'cornerLabelOverlay')) = [];
    if isempty(axs), error('Only overlay axes found in source figure.'); end
    % Choose the axes with most children (usually the primary plot)
    [~,i] = max(arrayfun(@(a) numel(allchild(a)), axs));
    ax = axs(i);
end

function cloneAxesProps(srcAx, dstAx)
    % Copy limits and labels text only (not positions)
    try, set(dstAx,'XLim',get(srcAx,'XLim'),'YLim',get(srcAx,'YLim')); end
    try, xlabel(dstAx, get(get(srcAx,'XLabel'),'String')); end
    try, ylabel(dstAx, get(get(srcAx,'YLabel'),'String')); end
    % Grid/scale/colormap if relevant
    try, grid(dstAx, gridState(srcAx)); end
    try, set(dstAx,'XScale',get(srcAx,'XScale'),'YScale',get(srcAx,'YScale')); end
    % Colormap only matters if images are present
    try, colormap(dstAx, colormap(srcAx)); end
end

function state = gridState(ax)
    gx = get(ax,'XGrid'); gy = get(ax,'YGrid');
    state = 'off';
    if (ischar(gx) && strcmpi(gx,'on')) || (ischar(gy) && strcmpi(gy,'on'))
        state = 'on';
    end
end
