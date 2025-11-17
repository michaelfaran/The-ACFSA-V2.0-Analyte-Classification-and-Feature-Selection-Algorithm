function figOut = collage1x3CopyPlus1Create(figs12, ...
    config, all_sensors_length, number_of_sensors, ...
    avg_dist_vec, Vektor_ARI, Gaussian_error_val,underlinetitle,lablesact,Lambda_1,Lambda_2)

% 1x3 collage, fixed canvas 17.8 cm × 4.45 cm (DPI-robust), zero gap:
% [ fig1 | fig2 | create3DAccuracyFigUP(...) ]
% - Copies axes(+legend/colorbar) preserving geometry in cm and color props.
% - Tile 3 is generated to the same cm-sized axes box as tile 2 (fallback inset).
% - White backgrounds forced. Panel labels A–C (top-left).
%
% Example:
% figOut = collage1x3CopyPlus1Create([figure(1) figure(2)], ...
%            config, all_sensors_length, 1:1:all_sensors_length, ...
%            avg_dist_vec, Vektor_ARI, Gaussian_error_val_mat);

    % ---------- inputs ----------
    assert(numel(figs12)==2 && all(arrayfun(@ishandle,figs12)), ...
        'Provide exactly two valid figure handles for tiles 1–2.');

    % ---------- fixed canvas ----------
    TOTAL_W = 17.7;   % cm
    TOTAL_H = 7/2*1.5;   % cm
    tileW   = TOTAL_W/3;
    tileH   = TOTAL_H;
    insetyy = 0.125;
    nnn=size(Gaussian_error_val,1);
    figOut = figure('Units','centimeters','Position',[2 2 TOTAL_W TOTAL_H], ...
        'Color','w','MenuBar','none','ToolBar','none', ...
        'NumberTitle','off','Name','Collage 1x3 (cm)', ...
        'PaperUnits','centimeters','PaperPositionMode','auto', ...
        'InvertHardcopy','off', ...
        'DefaultAxesUnits','centimeters', ...
        'DefaultAxesColor','w', ...
        'DefaultUipanelBackgroundColor','w', ...
        'DefaultUicontrolBackgroundColor','w');

    mkpanel = @(x0) uipanel('Parent',figOut,'Units','centimeters', ...
        'Position',[x0 0 tileW tileH], 'BorderType','none', ...
        'BackgroundColor','w','AutoResizeChildren','off');

    % Panels (1x3, zero gap)
    p1 = mkpanel(0);
    p2 = mkpanel(tileW);
    p3 = mkpanel(2*tileW);

    % ---------- copy figs -> panels ----------
    copyFigureIntoPanel_cm(figs12(1), p1,insetyy,nnn);
    copyFigureIntoPanel_cm(figs12(2), p2,insetyy,nnn);
% underlinedtitle('Yossi');  %underlinetitle
underlinedtitle(underlinetitle);
    % Reference axes box (cm) taken from tile 2, else inset fallback
     refAx  = findobj(p2,'Type','axes','-depth',1);
    % refAx=[];
    if ~isempty(refAx)
        refPosCM1 = getAxesPosCm(refAx(1));
        refPosCM2 = getAxesPosCm(refAx(2));
        refPosCM=refPosCM2;
        offset=0.06;
        gapY=0.3540-offset; %convereted from 0.12 in normaizlied here
        refPosCM(4)=refPosCM2(4)+refPosCM1(4)+1*gapY;
    else
        inset = 0.2;  % cm
        refPosCM = [inset inset tileW-2*inset tileH-2*inset];
    end

    % ---------- tile 3: create plot, match cm geometry ----------
    ax3 = axes('Parent',p3,'Units','centimeters','Position',refPosCM, ...
               'Color','w','Box','off');
    axes(ax3); %#ok<LAXES>

    try
        [ax1, ax2] = create3DAccuracyFigUP( ...
            config, all_sensors_length, number_of_sensors, ...
            avg_dist_vec, Vektor_ARI, Gaussian_error_val, ax3, p3, refPosCM,Lambda_1,Lambda_2);
    catch ME
        warning('create3DAccuracyFigUP failed: %s', ME.message);
        ax1 = ax3; ax2 = [];
    end

    % Force all axes in tile 3 to same cm box & white bg unless transparent
    axList3 = findobj(p3,'Type','axes');
    for a = reshape(axList3,1,[])
        set(a,'Units','centimeters','Position',refPosCM, ...
              'ActivePositionProperty','position','Box','off');
        col = get(a,'Color');
        if ~(ischar(col) && strcmp(col,'none')) && ...
           ~(isnumeric(col) && any(isnan(col)))
            set(a,'Color','w');
        end
    end

    % Legend in tile 3: top-center (horizontal)
    leg3 = findobj(p3,'Type','Legend');
    for L = reshape(leg3,1,[])
        % set(L,'Units','normalized','Position',[0.2 0.78 0.6 0.08], ...
        %       'Orientation','horizontal');
        hhh=get(findobj(p2,'Type','Legend'),'Position');
        hhh(1)=hhh(1)+0.05;
        hhh(3)= 0.57;
        if nnn==6
        hhh(1)=hhh(1)+0.05;
        elseif nnn==3
        hhh(1)=hhh(1)-0.125;           
        end
        % hhh(2)=hhh(2)-insety;
                set(L,'Units','normalized','Position',hhh, ...
              'Orientation','horizontal');
    end

    
    % ---------- panel labels ----------
    if lablesact==1
    placeCornerLabel_onPanel(p1,'A','topleft');
    placeCornerLabel_onPanel(p2,'B','topleft');
    placeCornerLabel_onPanel(p3,'C','topleft');
    end

    set(figOut,'PaperPositionMode','auto');  % export stays white

% ===== nested helpers =====
    function copyFigureIntoPanel_cm(srcFig, destPanel,insetyy,nnn)
    success = true;
    try
        axList = findobj(srcFig,'Type','axes','-depth',1);
        if isempty(axList), success=false; end
        Lall = findobj(srcFig,'Type','Legend','-depth',1);
        Call = findobj(srcFig,'Type','ColorBar','-depth',1);

        for ax = flipud(axList(:).')
            axPosCM = getAxesPosCm(ax);

            objs = ax;
            if ~isempty(Lall)
                m = find(arrayfun(@(h) isequal(get(h,'Axes'),ax), Lall),1);
                if ~isempty(m), objs(end+1)=Lall(m); end %#ok<AGROW>
            end
            if ~isempty(Call)
                m = find(arrayfun(@(h) isequal(get(h,'Axes'),ax), Call),1);
                if ~isempty(m), objs(end+1)=Call(m); end %#ok<AGROW>
            end

            newObjs = copyobj(objs, destPanel);

            % identify the axes among the copied objects
            isAx = arrayfun(@(h) strcmp(get(h,'Type'),'axes'), newObjs);
            if any(isAx)
                newAx = newObjs(find(isAx,1,'first'));
            else
                newAx = findobj(destPanel,'Type','axes','-depth',1); newAx = newAx(1);
            end

            % place geometry in cm and enforce white bg (unless transparent)
            set(newAx,'Units','centimeters','Position',axPosCM, ...
                      'ActivePositionProperty','position');
            try
                col = get(newAx,'Color');
                if ~(ischar(col) && strcmp(col,'none')) && ...
                   ~(isnumeric(col) && any(isnan(col)))
                    set(newAx,'Color','w');
                end
            end

            % clone color-related props
            cloneColorProps(ax, newAx);

            % relink colorbar to the new axes if present
            newCb = newObjs(arrayfun(@(h) strcmp(get(h,'Type'),'ColorBar'), newObjs));
            if ~isempty(newCb)
                try, set(newCb(1),'Axes',newAx); end %#ok<TRYNC>
            end
        end
    catch
        success = false;
    end

    if ~success
        % raster fallback at 1:1 panel size (handles yyaxis, etc.)
        tmp = [tempname,'.png'];
        exportgraphics(srcFig,tmp,'Resolution',150);
        img = imread(tmp); delete(tmp);
        axImg = axes('Parent',destPanel,'Units','centimeters', ...
                     'Position',[0 0 destPanel.Position(3) destPanel.Position(4)], ...
                     'Color','w','Visible','off');
        image(axImg,img); axis(axImg,'off','image');
    end
    % Lall.Position(2)=Lall.Position(2)-insetyy;
    set(findobj(gcf,'Type','Legend'),'Position',[  0.2332    0.7819-insetyy+0.02    0.5336    0.1261]);
    if nnn==3
    set(findobj(gcf,'Type','Legend'),'Position',[  0.2332+insetyy    0.7819-insetyy+0.02    0.3498    0.1261]);
    end

end

function poscm = getAxesPosCm(axh)
    old = get(axh,'Units');
    set(axh,'Units','centimeters');
    poscm = get(axh,'Position');
    set(axh,'Units',old);
end

function cloneColorProps(srcAx, dstAx)
    try, cmap = colormap(srcAx); if ~isempty(cmap), colormap(dstAx, cmap); end, end
    try, clim = get(srcAx,'CLim'); if ~isempty(clim), caxis(dstAx,clim); end, end
    try, alim = get(srcAx,'ALim'); if ~isempty(alim), set(dstAx,'ALim',alim); end, end
    try, cs   = get(srcAx,'ColorScale'); if ~isempty(cs), set(dstAx,'ColorScale',cs); end, end
end

function placeCornerLabel_onPanel(panelH, txt, where)
    % Overlay text so it stays on top regardless of children order.
    old = findobj(panelH,'Type','axes','Tag','cornerLabelOverlay');
    if ~isempty(old), delete(old); end
    ax = axes('Parent',panelH,'Units','normalized','Position',[0 0 1 1], ...
              'Color','none','XLim',[0 1],'YLim',[0 1], ...
              'Visible','off','HitTest','off','Tag','cornerLabelOverlay');
    uistack(ax,'top');
    insety = 0.125;
    inset=0.05;
    switch lower(where)
        case 'topleft',     x= inset; y=1-insety;  hal='left';  val='top';
        case 'topright',    x=1-inset; y=1-insety; hal='right'; val='top';
        case 'bottomleft',  x= inset; y= insety;   hal='left';  val='bottom';
        case 'bottomright', x=1-inset; y= insety;  hal='right'; val='bottom';
        otherwise, error('Unknown corner "%s"', where);
    end
    text(ax,x,y,txt,'Units','normalized','HorizontalAlignment',hal, ...
         'VerticalAlignment',val,'FontWeight','normal','FontSize',12, ...
         'Color','k','Interpreter','none');
end
end

function h = underlinedtitle(txt, fig)
    if nargin < 2 || isempty(fig), fig = gcf; end
    if iscell(txt), txt = txt{1}; end        % allow cellstr
    if isstring(txt), txt = txt(1); end      % string -> take first element
    txt = char(txt);                         % <-- remove quotes
    txt = escapeLatex(txt);                  % escape LaTeX specials

    % remove previous one on this figure
    delete(findobj(fig,'Type','textboxshape','Tag','underlinetitle'));
    % lll=fig.Position;
    % lll2=[lll(1)+lll(3)/2 lll(4)-2 2 1 ];
    % top-centered textbox, normalized to THIS figure
    % h = annotation(fig,'textbox',[0 0.96 1 0.04], ...
    %     'Units','normalized', ...
    %     'String', sprintf('\\underline{%s}', txt), ... % <- use sprintf
    %     'Interpreter','latex', ...
    %     'HorizontalAlignment','center','VerticalAlignment','top', ...
    %     'EdgeColor','none','Tag','underlinetitle');
t = txt;
if isstring(t), t = t(1); end
t = char(t);

sgtitle(sprintf('\\underline{%s}', escapeLatex(t)), 'Interpreter','latex','FontSize',10);


end

function s = escapeLatex(s)
    s = char(s);                          % ensure char out, not string
    s = strrep(s,'\','\textbackslash{}');
    s = strrep(s,'_','\_');  s = strrep(s,'%','\%');
    s = strrep(s,'$','\$');  s = strrep(s,'{','\{');
    s = strrep(s,'}','\}');  s = strrep(s,'#','\#');
    s = strrep(s,'&','\&');  s = strrep(s,'^','\^{}');
end

