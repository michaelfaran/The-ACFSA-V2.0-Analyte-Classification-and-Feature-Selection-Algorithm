function Main
% ACFSA V2 official GUI
% Two-step GUI:
%   Step 1  – Create/Open template and Import XLSX → produces DAT, analyte_name_vec, sensor_name_vec, n_measure
%   Step 2  – Configure ACFSA V2 options → calls:
%             Main_to_Call_GUI(supertitle, inputs_mat, analyte_name_vec, sensor_name_vec, DAT, n_measure)
%
% inputs_mat = [inflateFlag, artificialSTDmultiplier, classifierFlag, weightedFSFlag, allowedErrWP1_pct, allowedErrWP2_pct]
%   - inflateFlag:               0/1
%   - artificialSTDmultiplier:   0 = no artificial dataset; >0 = STD×multiplier stress-test (must be ≥ 0)
%   - classifierFlag:            0 = QDA, 1 = Voronoi
%   - weightedFSFlag:            0/1
%   - allowedErrWP1_pct:         working point 1, allowed classification error (%) per sensor in minimal dataset (default 3)
%   - allowedErrWP2_pct:         working point 2, allowed classification error (%) per sensor in minimal dataset (default 0.5)

    % ---- Local state (for Step 2) ----
    DAT_local                 = [];
    analyte_name_vec_local    = {};
    sensor_name_vec_local     = {};
    n_measure_local           = [];
    canProceed                = false;

    % ---- STEP 1 UI ----
    fig = uifigure('Name','Measurement XLSX: Template & Import', ...
                   'Position',[80 80 1000 640]); % enlarged window
    gl  = uigridlayout(fig,[6 1]);
    gl.RowHeight = {110, 60, '2x', 60, 60, 44};

    % Inputs
    pin = uipanel(gl,'Title','Template settings','FontWeight','bold');
    g = uigridlayout(pin,[2 6]);
    g.RowHeight   = {28, '1x'};
    g.ColumnWidth = {150, 100, 150, 100, 170, '1x'};

    uilabel(g,'Text','# Analytes:','HorizontalAlignment','right');
    edtA = uieditfield(g,'numeric','Limits',[1 Inf],'RoundFractionalValues','on','Value',3);

    uilabel(g,'Text','# Sensors:','HorizontalAlignment','right');
    edtS = uieditfield(g,'numeric','Limits',[1 Inf],'RoundFractionalValues','on','Value',20);

    uilabel(g,'Text','# Measurements / analyte:','HorizontalAlignment','right');
    edtM = uieditfield(g,'numeric','Limits',[1 Inf],'RoundFractionalValues','on','Value',3);

    % Path (template create/open target)
    pPath = uipanel(gl,'Title','Template file path (for Create/Open)','FontWeight','bold');
    pathLabel = uilabel(pPath,'Text',defaultTemplatePath(), ...
        'Interpreter','none','Position',[10 5 pPath.InnerPosition(3)-20 25]);
    addlistener(pPath,'SizeChanged',@(~,~)set(pathLabel,'Position',[10 5 pPath.InnerPosition(3)-20 25]));

    % Instructions
    pHelp = uipanel(gl,'Title','Instructions','FontWeight','bold');
    hgl  = uigridlayout(pHelp,[1 1]); hgl.RowHeight = {'1x'}; hgl.ColumnWidth = {'1x'};
    helpLbl = uilabel(hgl,'WordWrap','on','HorizontalAlignment','left','FontSize',12);

    % Last imported file path
    pImp = uipanel(gl,'Title','Last chosen import file','FontWeight','bold');
    importPathLabel = uilabel(pImp,'Text','(none chosen yet)', ...
        'Interpreter','none','Position',[10 5 pImp.InnerPosition(3)-20 25]);
    addlistener(pImp,'SizeChanged',@(~,~)set(importPathLabel,'Position',[10 5 pImp.InnerPosition(3)-20 25]));

    % Buttons
    row = uigridlayout(gl,[1 3]); row.ColumnWidth = {260, 260, 260};
    btnCreateOpen = uibutton(row,'Text','Create Template & Open in Excel', ...
        'FontWeight','bold','ButtonPushedFcn',@onCreateAndOpen);
    btnImport     = uibutton(row,'Text','Import from XLSX (choose file)', ...
        'FontWeight','bold','ButtonPushedFcn',@onImport);
    btnNext       = uibutton(row,'Text','Next: Configure & Run ACFSA V2', ...
        'FontWeight','bold','Enable','off','ButtonPushedFcn',@onProceed);

    % Status
    pStatus = uipanel(gl,'Title','Status','FontWeight','bold');
    sgl = uigridlayout(pStatus,[1 1]); sgl.RowHeight = {'1x'}; sgl.ColumnWidth = {'1x'};
    status  = uilabel(sgl,'Text','Ready','FontColor',[0.1 0.45 0.1], 'HorizontalAlignment','left');

    % ---- Hook instruction auto-update ----
    edtS.ValueChangedFcn = @updateInstr;
    edtM.ValueChangedFcn = @updateInstr;
    updateInstr();

    % ---- Callbacks (STEP 1) ----
    function onCreateAndOpen(~,~)
        A = chkInt(edtA.Value,'# Analytes',1); if isnan(A), return; end
        S = chkInt(edtS.Value,'# Sensors',1); if isnan(S), return; end
        M = chkInt(edtM.Value,'# Measurements / analyte',1); if isnan(M), return; end

        file = defaultTemplatePath();
        try
            header = [{'Label'}, arrayfun(@(k)sprintf('sensor %d',k), 1:S, 'UniformOutput', false)];
            labels = buildLabels(A,M);                % AxM rows
            body   = num2cell(nan(A*M, S));           % numeric placeholders
            out    = [header; [labels body]];

            if exist(file,'file')
                delete(file); % overwrite cleanly
            end
            writecell(out, file, 'FileType','spreadsheet');

            setStatus(sprintf('Template written: %s', file), true);
            tryOpen(file);
            pathLabel.Text = file;
        catch ME
            uialert(fig,ME.message,'Write/Open Error','Icon','error');
            setStatus('Error creating or opening file.', false);
        end
    end

    function onImport(~,~)
        [fname, fpath] = uigetfile({'*.xlsx;*.xls','Excel files (*.xlsx, *.xls)'}, ...
                                   'Choose measurement XLSX to import', pwd);
        if isequal(fname,0)
            setStatus('Import canceled by user.', false);
            return;
        end
        file = fullfile(fpath, fname);
        importPathLabel.Text = file;

        try
            T = readtable(file,'TextType','string','PreserveVariableNames',true);
        catch ME
            uialert(fig,ME.message,'Read Error','Icon','error');
            setStatus('Read error.', false); return;
        end
        if width(T) < 2
            uialert(fig,'Expected at least 2 columns: a Label column + ≥1 sensor columns.','Format Error','Icon','warning');
            setStatus('Bad format.', false); return;
        end

        info = detectFromTable(T);
        if ~info.ok
            uialert(fig,info.errMsg,'Format Error','Icon','error');
            setStatus('Bad format.', false); return;
        end

        msg = composeSummaryMessage(info);
        choice = uiconfirm(fig, msg, 'Review detected structure', ...
            'Options',{'Use detected','Enter manually','Cancel'}, ...
            'DefaultOption',1,'CancelOption',3, 'Icon','info');

        switch choice
            case 'Cancel'
                setStatus('Import canceled by user.', false);
                return;

            case 'Enter manually'
                prompt = {'# Analytes (A):', '# Sensors (S):', '# Measurements per analyte (M):'};
                defans = {num2str(numel(info.analyte_name_vec)), num2str(numel(info.sensor_name_vec)), num2str(info.n_measure)};
                answer = inputdlg(prompt,'Override expected sizes', [1 45], defans);
                if isempty(answer)
                    setStatus('Import canceled by user.', false);
                    return;
                end
                Aexp = str2double(answer{1});
                Sexp = str2double(answer{2});
                Mexp = str2double(answer{3});
                if any(~isfinite([Aexp Sexp Mexp])) || any([Aexp Sexp Mexp] < 1) || ...
                   any(abs([Aexp Sexp Mexp]-round([Aexp Sexp Mexp]))>eps)
                    uialert(fig,'All expected sizes must be positive integers.','Input Error','Icon','warning');
                    setStatus('Bad expected sizes.', false); return;
                end

                analyte_name_vec = arrayfun(@(k)sprintf('analyte_name_%d',k), 1:Aexp, 'UniformOutput', false);
                sensor_name_vec  = arrayfun(@(k)sprintf('sensor_name_%d',k),  1:Sexp, 'UniformOutput', false);
                n_measure        = round(Mexp);
                DAT              = zeros(Aexp * n_measure, Sexp);

                assignin('base','analyte_name_vec',analyte_name_vec(:).');
                assignin('base','sensor_name_vec', sensor_name_vec(:).');
                assignin('base','n_measure',       n_measure);
                assignin('base','DAT',             DAT);

                analyte_name_vec_local = analyte_name_vec(:).';
                sensor_name_vec_local  = sensor_name_vec(:).';
                n_measure_local        = n_measure;
                DAT_local              = DAT;
                canProceed             = true;
                btnNext.Enable         = 'on';

                uialert(fig, sprintf(['Placeholders created in the MATLAB workspace:\n' ...
                    ' • analyte_name_vec: 1×%d cell\n' ...
                    ' • sensor_name_vec : 1×%d cell\n' ...
                    ' • n_measure       : %d\n' ...
                    ' • DAT             : %d×%d double\n\n' ...
                    'Please edit their values in the MATLAB workspace. Then, you can now click "Next" to configure ACFSA V2.'], ...
                    Aexp, Sexp, n_measure, size(DAT,1), size(DAT,2)), ...
                    'Placeholders created','Icon','info');

                setStatus(sprintf('Manual placeholders created: rows=%d, sensors=%d, analytes=%d, n_measure=%d', ...
                    size(DAT,1), Sexp, Aexp, n_measure), true);
                return;

            otherwise % 'Use detected'
                Aexp = numel(info.analyte_name_vec);
                Sexp = numel(info.sensor_name_vec);
                Mexp = info.n_measure;

                [ok, softWarn, hardErr] = verifyAgainstExpectations(T, info, Aexp, Sexp, Mexp);
                if ~ok && ~isempty(hardErr)
                    uialert(fig, hardErr, 'Format Error','Icon','error');
                    setStatus('Failed verification.', false); return;
                end
                if ~isempty(softWarn)
                    proceed = uiconfirm(fig, softWarn + newline + "Continue import?", ...
                        'Warnings', 'Options',{'Continue','Cancel'}, 'DefaultOption',1, 'CancelOption',2, 'Icon','warning');
                    if strcmp(proceed,'Cancel')
                        setStatus('Import canceled due to warnings.', false);
                        return;
                    end
                end

                NS = numel(info.sensor_name_vec);
                NR = height(T);
                DAT = nan(NR, NS);
                for c = 1:NS
                    col = T.(c+1);
                    if isnumeric(col)
                        DAT(:,c) = double(col);
                    else
                        DAT(:,c) = str2double(string(col)); %#ok<ST2NM>
                    end
                end

                analyte_name_vec = info.analyte_name_vec(:).';
                sensor_name_vec  = info.sensor_name_vec(:).';
                n_measure        = info.n_measure;

                assignin('base','DAT',DAT);
                assignin('base','analyte_name_vec',analyte_name_vec);
                assignin('base','sensor_name_vec',sensor_name_vec);
                assignin('base','n_measure',n_measure);

                DAT_local              = DAT;
                analyte_name_vec_local = analyte_name_vec;
                sensor_name_vec_local  = sensor_name_vec;
                n_measure_local        = n_measure;
                canProceed             = true;
                btnNext.Enable         = 'on';

                setStatus(sprintf('Imported OK: rows=%d, sensors=%d, analytes=%d, n_measure=%d', ...
                    size(DAT,1), numel(sensor_name_vec), numel(analyte_name_vec), n_measure), true);
        end
    end

    function onProceed(~,~)
        if ~canProceed || isempty(DAT_local) || isempty(analyte_name_vec_local) ...
                || isempty(sensor_name_vec_local) || isempty(n_measure_local)
            uialert(fig,'Please import (or create) data first.','Missing Data','Icon','warning');
            return;
        end

        % ---- STEP 2 UI ----
        fig2 = uifigure('Name','ACFSA V2: Configure & Run','Position',[120 120 800 520]);

        % Layout for added working-point inputs + info text (no checkbox)
        gl2 = uigridlayout(fig2,[9 2]);
        gl2.RowHeight   = {38, 38, 38, 38, 38, 38, 38, 40, '1x'};
        gl2.ColumnWidth = {360, '1x'};

        % Title
        uilabel(gl2,'Text','Dataset title (supertitle):','HorizontalAlignment','right');
        edtTitle = uieditfield(gl2,'text','Value','ACFSA V2 Run');

        % Inflate uncertainty? (0/1)
        uilabel(gl2,'Text','Inflate data to cover statistical uncertainty? (0/1):','HorizontalAlignment','right');
        edtInflate = uieditfield(gl2,'numeric','Limits',[0 1],'Value',1,'RoundFractionalValues','off');

        % Artificial dataset (numeric multiplier, 0 = none)
        uilabel(gl2,'Text','Artificial dataset: STD multiplier (0 = none, ≥ 0):','HorizontalAlignment','right');
        edtArtificial = uieditfield(gl2,'numeric','Limits',[0 Inf],'Value',0, ...
            'Tooltip','0 disables artificial dataset. >0 multiplies class STD to synthesize a stress-test set.');

        % Classifier: 0 = QDA, 1 = Voronoi
        uilabel(gl2,'Text','Classifier (0 = QDA, 1 = Voronoi):','HorizontalAlignment','right');
        edtClassifier = uieditfield(gl2,'numeric','Limits',[0 1],'Value',0,'RoundFractionalValues','off');

        % Weighted feature selection? (0/1)
        uilabel(gl2,'Text','Weighted feature selection? (0 = weighted, 1 = standard):','HorizontalAlignment','right');
        edtWeighted = uieditfield(gl2,'numeric','Limits',[0 1],'Value',1,'RoundFractionalValues','off');

        % Working points – labels updated to say "per sensor"
        uilabel(gl2,'Text','Allowed classification error per sensor (%) — Working point 1:','HorizontalAlignment','right');
        edtWP1 = uieditfield(gl2,'numeric','Limits',[0 100],'Value',0.5, ...
            'Tooltip','Working point 1: allowed classification error per sensor (%) in the minimal dataset.');

        uilabel(gl2,'Text','Allowed classification error per sensor (%) — Working point 2:','HorizontalAlignment','right');
        edtWP2 = uieditfield(gl2,'numeric','Limits',[0 100],'Value',3, ...
            'Tooltip','Working point 2 (default): allowed classification error per sensor (%) in the minimal dataset.');

        % Informative note (no checkbox)
        infoLbl = uilabel(gl2,'Text', ...
            'This GUI will be closed before running. Input parameters will be printed again at the end of the run.');
        infoLbl.Layout.Row = 8; infoLbl.Layout.Column = [1 2];

        % Run button
        runRow = uigridlayout(gl2,[1 2]); runRow.Layout.Row = 9; runRow.Layout.Column = [1 2];
        btnRun = uibutton(runRow,'Text','Run ACFSA','FontWeight','bold','ButtonPushedFcn',@onRunACFSA);
        uilabel(runRow,'Text',''); % filler

        function onRunACFSA(~,~)
            supertitle = edtTitle.Value;

            % Validate inputs
            aSTD = edtArtificial.Value;
            if ~(isfinite(aSTD) && aSTD >= 0)
                uialert(fig2,'Artificial STD multiplier must be a finite number ≥ 0.','Input Error','Icon','warning');
                return;
            end

            infl = edtInflate.Value;
            clsf = edtClassifier.Value;
            wfs  = edtWeighted.Value;
            wp1  = edtWP1.Value;
            wp2  = edtWP2.Value;

            if ~(isfinite(infl) && any(infl == [0 1]))
                uialert(fig2,'Inflate flag must be exactly 0 or 1.','Input Error','Icon','warning'); return;
            end
            if ~(isfinite(clsf) && any(clsf == [0 1]))
                uialert(fig2,'Classifier must be 0 (QDA) or 1 (Voronoi).','Input Error','Icon','warning'); return;
            end
            if ~(isfinite(wfs) && any(wfs == [0 1]))
                uialert(fig2,'Weighted feature selection flag must be exactly 0 or 1.','Input Error','Icon','warning'); return;
            end
            if ~(isfinite(wp1) && wp1 >= 0 && wp1 <= 100)
                uialert(fig2,'Working point 1 must be a percentage in [0, 100].','Input Error','Icon','warning'); return;
            end
            if ~(isfinite(wp2) && wp2 >= 0 && wp2 <= 100)
                uialert(fig2,'Working point 2 must be a percentage in [0, 100].','Input Error','Icon','warning'); return;
            end

            % Anchor classic figure (legacy plotting friendliness)
            anchorFig = figure('Visible','off','HandleVisibility','on','IntegerHandle','off', ...
                               'NumberTitle','off','Name','ACFSA-AnchorClassic');
            try, set(0,'CurrentFigure',anchorFig); catch, end
            drawnow;

            % ALWAYS close the GUIs before running (no checkbox)
            try, if isvalid(fig2), delete(fig2); end, catch, end
            try, if isvalid(fig),  delete(fig);  end, catch, end

            % Build inputs and run (6-D)
            inputs_mat = double([infl, aSTD, clsf, wfs, wp1, wp2]); %#ok<NASGU>
            input_vec  = inputs_mat; %#ok<NASGU>

            try
                Main_to_Call_GUI(supertitle, inputs_mat, ...
                    analyte_name_vec_local, sensor_name_vec_local, DAT_local, n_measure_local);

                % Print parameters to Command Window at end of run
                fprintf('\n[ACFSA V2] Run completed.\n');
                fprintf('Title: %s\n', supertitle);
                fprintf('Inflate uncertainty: %g\n', infl);
                fprintf('Artificial STD multiplier: %g\n', aSTD);
                fprintf('Classifier: %s\n', tern(clsf==0,'QDA','Voronoi'));
                fprintf('Weighted feature selection: %s\n', tern(wfs==1,'ON','OFF'));
                fprintf('Allowed classification error per sensor (%%) — WP1: %.4g, WP2: %.4g\n', wp1, wp2);
                fprintf('Dataset: %d rows × %d sensors; analytes: %d; n_measure: %d\n', ...
                    size(DAT_local,1), size(DAT_local,2), numel(analyte_name_vec_local), n_measure_local);
            catch ME
                % Even on error, print the chosen parameters for traceability
                fprintf('\n[ACFSA V2] Run ERROR: %s\n', ME.message);
                fprintf('Parameters were: Inflate=%g, aSTD=%g, Classifier=%s, WeightedFS=%s, WP1=%.4g%%, WP2=%.4g%%\n', ...
                    infl, aSTD, tern(clsf==0,'QDA','Voronoi'), tern(wfs==1,'ON','OFF'), wp1, wp2);
                try, errordlg(ME.message,'Run Error'); end
            end
        end
    end

    % ---- Helpers ----
    function updateInstr(~,~)
        S = round(edtS.Value);
        M = round(edtM.Value);
        helpLbl.Text = sprintf([ ...
            'Please complete the Excel template before importing:\n' ...
            '• Add sensor names in row 1 and analyte names in the first column.\n' ...
            '• Enter the numeric measurement values in the data area (no text in sensor columns).\n' ...
            '• For each analyte, its measurements must occupy consecutive rows (exactly %d rows per analyte).\n' ...
            '• Each measurement row must contain one numeric value per sensor (row length = %d sensor columns).' ...
            '\n Notes:\n'...
            '• This algorithm assumes that the measurements are already normalized sensor responses.\n'...
            '• The classification-error output and any synthetic dataset (if generated) assume that each analyte’s measurements follow a Gaussian distribution.'...
            ], M, S);
    end

    function setStatus(txt, ok)
        status.Text = txt;
        status.FontColor = iff(ok, [0.1 0.45 0.1], [0.6 0 0]);
    end

    function p = defaultTemplatePath()
        p = fullfile(pwd,'measurement_data.xlsx');
    end

    function labels = buildLabels(A,M)
        R = A*M; labels = cell(R,1); k = 0;
        for a = 1:A
            for m = 1:M
                k = k + 1;
                labels{k} = sprintf('analyte %d, relative fluorescence peak measurement %d', a, m);
            end
        end
    end

    function n = chkInt(v, name, minv)
        if isempty(v) || ~isfinite(v) || v < minv || abs(v-round(v)) > eps
            uialert(fig, sprintf('%s must be an integer ≥ %d.', name, minv), 'Invalid Input','Icon','warning');
            n = NaN;
        else
            n = round(v);
        end
    end

    function tryOpen(fname)
        if ispc, winopen(fname); return; end
        if ismac, system(['open "', fname, '" &']); return; end
        system(['xdg-open "', fname, '" &']);
    end

    function y = iff(c,a,b), if c, y=a; else, y=b; end, end
    function y = tern(cond,a,b), if cond, y=a; else, y=b; end, end

    % ---- Detection / Verification utilities ----
    function info = detectFromTable(T)
        info = struct('ok',true,'errMsg',"", ...
            'labelVar',"",'labelCol',[], ...
            'analyte_name_vec',{{}}, 'counts',[], 'n_measure',[], ...
            'sensor_name_vec',{{}}, 'NS',[], 'NR',[], 'DAT_raw',[]);
        try
            info.NR = height(T);
            info.labelVar = T.Properties.VariableNames{1};
            info.labelCol = string(T.(info.labelVar));

            analyteAll = strip(extractBefore(info.labelCol, ","));
            noCommaIdx = analyteAll == "";
            analyteAll(noCommaIdx) = strip(info.labelCol(noCommaIdx));

            [~, firstIdx] = unique(analyteAll,'stable');
            info.analyte_name_vec = cellstr(analyteAll(sort(firstIdx)).');  % 1xA

            [grp,~] = grp2idx(analyteAll);
            info.counts = accumarray(grp,1);
            info.n_measure = mode(info.counts);

            info.sensor_name_vec = cellstr(T.Properties.VariableNames(2:end));
            info.NS = numel(info.sensor_name_vec);

            info.DAT_raw = nan(info.NR, info.NS);
            for c = 1:info.NS
                col = T.(c+1);
                if isnumeric(col)
                    info.DAT_raw(:,c) = double(col);
                else
                    info.DAT_raw(:,c) = str2double(string(col));
                end
            end

            if info.NS < 1
                info.ok = false; info.errMsg = "Expected at least one sensor column."; return;
            end
            if any(strlength(string(info.sensor_name_vec))==0)
                info.ok = false; info.errMsg = "Empty sensor header detected."; return;
            end
            if numel(unique(lower(string(info.sensor_name_vec)))) < info.NS
                info.ok = false; info.errMsg = "Duplicate sensor headers detected (case-insensitive)."; return;
            end
        catch ME
            info.ok = false; info.errMsg = ME.message;
        end
    end

    function msg = composeSummaryMessage(info)
        A = numel(info.analyte_name_vec);
        S = numel(info.sensor_name_vec);
        M = info.n_measure;
        countsTxt = sprintf('%s', strjoin(string(info.counts.'),', '));
        msg = sprintf([ ...
            'Detected from file:\n' ...
            ' • Analytes (A): %d  →  %s\n' ...
            ' • Sensors  (S): %d  →  %s\n' ...
            ' • Measurements per analyte (mode): %d\n' ...
            ' • Row counts per analyte: [%s]\n\n' ...
            'Use these, or enter manually?'], ...
            A, strjoin(info.analyte_name_vec, ', '), ...
            S, strjoin(info.sensor_name_vec, ', '), ...
            M, countsTxt);
    end

    function [ok, softWarn, hardErr] = verifyAgainstExpectations(T, info, Aexp, Sexp, Mexp)
        ok = true; softWarn = ""; hardErr = "";

        if width(T)-1 ~= Sexp
            ok = false;
            hardErr = sprintf('Header sensor count (%d) does not match expected S (%d).', width(T)-1, Sexp);
            return;
        end
        if numel(info.analyte_name_vec) ~= Aexp
            ok = false;
            hardErr = sprintf('Detected analyte classes (%d) do not match expected A (%d).', numel(info.analyte_name_vec), Aexp);
            return;
        end
        if any(all(isnan(info.DAT_raw),1))
            badCols = find(all(isnan(info.DAT_raw),1));
            ok = false;
            hardErr = sprintf('One or more sensor columns are entirely non-numeric/NaN. Columns: %s', mat2str(badCols));
            return;
        end

        if ~all(info.counts==Mexp)
            softWarn = softWarn + sprintf('Non-uniform or unexpected row counts per analyte. Expected M=%d, detected: [%s].\n', ...
                                          Mexp, strjoin(string(info.counts.'),', '));
        end

        % --- Label-format check (unchanged) ---
        lbl = string(info.labelCol);
        lbl = regexprep(lbl, '[\x00-\x1F\x7F\xA0\u2000-\u200B\u2028\u2029\u202F\u205F\u3000]', ' ');
        lbl = regexprep(lbl, '\s+', ' ');
        lbl = strtrim(lbl);
        measureRow = any(isfinite(info.DAT_raw), 2);
        idxData    = find(measureRow);
        tf = false(size(lbl));
        if ~isempty(idxData)
            L  = cellstr(lbl(idxData));
            mskHasWord = ~cellfun('isempty', regexp(L, 'measur(e)?ment', 'once', 'ignorecase'));
            patt = '^[^,]+,\s*.*measur(e)?ment\s*[0-9]+\s*$';
            subIdx = idxData(mskHasWord);
            Lsub   = L(mskHasWord);
            tf(subIdx) = ~cellfun('isempty', regexp(Lsub, patt, 'once', 'ignorecase'));
            if ~isempty(subIdx)
                ratio = mean(tf(subIdx));
                if ratio < 0.6
                    bad = subIdx(~tf(subIdx));
                    preview = string(lbl(bad(1:min(3,numel(bad)))));
                    softWarn = softWarn + ...
                        "Label format diverges from template ('<name>, measurement j'). Examples: " + ...
                        strjoin(preview, " | ") + newline;
                end
            end
        end

        if any(startsWith(string(T.Properties.VariableNames(2:end)),"Var"))
            softWarn = softWarn + "Sensor headers look like default ''VarX''. Consider renaming columns.\n";
        end
        if numel(unique(string(T.Properties.VariableNames(2:end)))) < (width(T)-1)
            softWarn = softWarn + "Sensor headers duplicated with case-sensitive comparison.\n";
        end
    end
end
