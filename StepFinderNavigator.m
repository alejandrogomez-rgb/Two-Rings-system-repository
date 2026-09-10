function StepFinderNavigator()
%STEPFINDERNAVIGATOR Interactive browser + detector for conductance steps
% (plateaus) in fast-breaking MCBJ traces.
%
% This is built on the same navigation scheme as FastBreakingNavigator.m:
% it steps through the breaking segment of a sequence of "..._N.dat"
% fast-breaking files one at a time, but on top of that it automatically
% detects conductance "steps" (plateaus) in each trace and overlays them
% so you can quickly eyeball whether the detection is doing a good job.
% Once you're happy with the detection (tune the parameters below and
% re-run with F5 if not), press 'e' to export .txt files, grouped by how
% many steps each trace has.
%
% ---- How a step is detected -------------------------------------------
%   1. d [nm] and y = log10(G/G0) are computed exactly as in
%      loadRawData.m / FastBreakingNavigator.m.
%   2. A Savitzky-Golay filter (dependency-free, no toolbox needed) is
%      used to get a smoothed local derivative dy/dd [decades of G0 / nm].
%   3. Points where |dy/dd| is small (below Parameters.DerivThresh) AND
%      y corresponds to Parameters.StepGmin < G/G0 < Parameters.StepGmax
%      are flagged as "quiet" (flat).
%   4. Consecutive quiet points are grouped into runs; short interruptions
%      (noise blips that briefly push the derivative over threshold) are
%      bridged if the gap is below Parameters.MergeGapNm.
%   5. A run is kept as a genuine "step" only if its displacement extent
%      falls within [Parameters.StepLenMin, Parameters.StepLenMax] nm.
%
% USAGE:
%   StepFinderNavigator()            % prompts you to pick a .dat file
%   StepFinderNavigator(fullpath)    % start directly from a known file
%
% Any file matching the same "<prefix>_<number>.dat" naming as the file
% you pick/pass will be included, sorted by that trailing number.
%
% CONTROLS (click on the navigator figure first so it has focus):
%   Right arrow / d   -> next trace
%   Left arrow  / a   -> previous trace
%   Spacebar          -> toggle EXCLUDE current trace from export
%   c                 -> clear all exclusions
%   g                 -> jump to a specific trace number
%   e                 -> export grouped .txt files for ALL traces
%   q                 -> close the navigator
%
% ---------------------------------------------------------------
% ---- EDIT THESE IF YOUR SETUP CHANGES (same as loadRawData.m) ----
Parameters.Vtodconversion = 0.3e-6;     % m/V, piezo voltage-to-displacement conversion
Parameters.attenuation    = 8e-5;       % dimensionless reduction/attenuation factor
Parameters.ZeroDisp       = 0.3;        % G/G0 threshold used to zero-reference displacement
Parameters.minlvl         = 1e-8;       % conductance floor before log10 (no baseline offset)
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---- STEP-DETECTION PARAMETERS: tune these while browsing -----
Parameters.SGWindowNm   = 0.025;   % nm, total width of the SG derivative window
Parameters.SGPolyOrder  = 2;      % polynomial order for the SG filter
Parameters.DerivThresh  = 5;     % decades of G0 per nm; below this counts as "flat"
Parameters.MergeGapNm   = 0.025;   % nm; bridges brief noise-driven interruptions
Parameters.StepGmin     = 5e-6;   % G0; steps are only ever considered above this
Parameters.StepGmax     = 5e-1;   % G0; steps are only ever considered below this
Parameters.StepLenMin   = 0.05;   % nm; minimum accepted step length
Parameters.StepLenMax   = 0.5;   % nm; maximum accepted step length
%   -> Too few / fragmented steps found?  Increase DerivThresh or MergeGapNm.
%   -> Steps blurring together / noise picked up as steps? Decrease them,
%      or shrink SGWindowNm so edges are less smeared.
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---- EXPORT PARAMETERS -----------------------------------------
Parameters.ExportMinSteps = 3;    % smallest step-count group to export
Parameters.ExportMaxSteps = 8;    % largest step-count group to export
% Only step counts in [ExportMinSteps, ExportMaxSteps] that have at
% least one matching (non-excluded) trace will produce a .txt file.
% ---------------------------------------------------------------

    % ---- 1) Pick a representative file, if not given ----
    if nargin < 1 || isempty(startFile)
        [fname, fpath] = uigetfile('*.dat', 'Select any .dat file from the fast-breaking series');
        if isequal(fname, 0)
            disp('No file selected. Exiting.');
            return
        end
        startFile = fullfile(fpath, fname);
    end
    if ~exist(startFile, 'file')
        error('File not found: %s', startFile);
    end

    % ---- 2) Derive the "<prefix>_" pattern and gather+sort matching files ----
    [folderPath, baseName, ~] = fileparts(startFile);
    parts = regexp(baseName, '_', 'split');
    if numel(parts) > 1
        prefix = strjoin(parts(1:end-1), '_');
    else
        prefix = baseName;
    end

    dList = dir(fullfile(folderPath, [prefix, '_*.dat']));
    if isempty(dList)
        error('No files found matching pattern "%s_*.dat" in %s', prefix, folderPath);
    end

    numericIDs = zeros(numel(dList), 1);
    for i = 1:numel(dList)
        nameNoExt = dList(i).name(1:end-4);
        p2 = regexp(nameNoExt, '_', 'split');
        val = str2double(p2{end});
        if isnan(val)
            val = Inf;  % non-numeric suffixes go last
        end
        numericIDs(i) = val;
    end
    [~, sortIdx] = sort(numericIDs);
    dList = dList(sortIdx);

    Nfiles = numel(dList);
    fileList = fullfile({dList.folder}, {dList.name});
    fprintf('Found %d fast-breaking trace file(s) matching "%s_*.dat"\n', Nfiles, prefix);

    % ---- 3) Read the breaking speed once, from the first file ----
    speed_const = readBreakingSpeed(fileList{1}, Parameters);

    % ---- 4) Build the navigator figure ----
    fig = figure('Color', 'w', 'Name', 'Step Finder Navigator', ...
                 'NumberTitle', 'off', 'KeyPressFcn', @keyPress);
    ax = axes(fig);
    set(ax, 'Color', 'w');

    data.fileList     = fileList;
    data.Nfiles        = Nfiles;
    data.currentIndex  = 1;
    data.excluded      = false(1, Nfiles);
    data.cache         = cell(Nfiles, 1);   % each entry: [d, y, g]
    data.stepCache     = cell(Nfiles, 1);   % each entry: steps struct array
    data.speed_const   = speed_const;
    data.Parameters    = Parameters;
    data.ax            = ax;
    data.folderPath    = folderPath;
    data.prefix        = prefix;
    guidata(fig, data);

    fprintf(['\nControls: Right/d = next, Left/a = previous, Space = exclude/include, ' ...
             'c = clear exclusions, g = go to #, e = export grouped .txt, q = quit\n\n']);

    renderCurrent(fig);
end

%% ========================================================================
function speed_const = readBreakingSpeed(fullFilename, Parameters)
%READBREAKINGSPEED Extract "Breaking speed 2" (V/s) from the file header
% and convert it to a displacement-rate constant in m/s.

    fid = fopen(fullFilename, 'r');
    fc  = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    fc  = fc{:};

    startline = find(strncmp('@', fc, 1), 1, 'first');
    if isempty(startline)
        header = fc;
    else
        header = fc(1:startline-1);
    end

    idx = find(strncmp('Breaking speed 2', header, 16));
    if isempty(idx)
        error('Could not find "Breaking speed 2" in the header of %s', fullFilename);
    end
    piezo_speed  = sscanf(header{idx(1)}, 'Breaking speed 2 : %f');  % V/s
    speed_const  = piezo_speed * Parameters.Vtodconversion;          % m/s
end

%% ========================================================================
function [d, y, g, ok] = readBreakingTrace(fullFilename, speed_const, Parameters)
%READBREAKINGTRACE Read the breaking segment (between the 1st and 2nd "@")
% of one .dat file and return displacement d [nm], log10(G/G0) y, and the
% (floor-clipped) conductance g [G0] used to compute y.

    d = []; y = []; g = []; ok = false;

    fid = fopen(fullFilename, 'r');
    fc  = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    fc  = fc{:};

    atLines = find(strncmp('@', fc, 1));
    if numel(atLines) < 2
        warning('%s: fewer than 2 "@" sections, skipping.', fullFilename);
        return
    end

    % Detect number of columns from the first data line
    tmp  = regexp(fc(atLines(1)+1), '\s*([\S]*)\s*', 'tokens');
    cols = numel(tmp{1});

    fid = fopen(fullFilename, 'r');
    if cols == 2
        data = textscan(fid, '%f %f', 'HeaderLines', atLines(1));
    else
        data = textscan(fid, '%f %f %f', 'HeaderLines', atLines(1));
    end
    fclose(fid);

    t = data{1};
    g = data{2};

    if isempty(t)
        warning('%s: empty breaking segment, skipping.', fullFilename);
        return
    end

    d = t * speed_const * Parameters.attenuation * 1e9;   % nm

    % Zero-reference displacement at first point where g < ZeroDisp
    try
        idx0   = find(g < Parameters.ZeroDisp);
        idx0   = idx0(1);
        offset = d(idx0);
        d      = d - offset;
    catch
        % trace never drops below ZeroDisp -> leave d unshifted
    end

    g(g < Parameters.minlvl) = Parameters.minlvl;
    y  = log10(g);
    ok = true;
end

%% ========================================================================
function C = sgolayCoeffs(halfWin, polyOrder, dspacing)
%SGOLAYCOEFFS Dependency-free Savitzky-Golay coefficient matrix.
% Row 1 gives smoothing coefficients, row 2 gives 1st-derivative
% coefficients (in units per unit of dspacing), etc.

    idx = (-halfWin:halfWin)';
    x   = idx * dspacing;
    A   = ones(numel(x), polyOrder + 1);
    for p = 1:polyOrder
        A(:, p+1) = x.^p;
    end
    C = (A' * A) \ A';
end

%% ========================================================================
function dy = sgDerivative(y, dspacing, halfWin, polyOrder)
%SGDERIVATIVE Smoothed first derivative dy/dspacing via Savitzky-Golay.
% Edge points (within halfWin of either end) are returned as NaN.

    n  = numel(y);
    dy = nan(n, 1);
    if n < (2*halfWin + 1)
        return
    end
    C = sgolayCoeffs(halfWin, polyOrder, dspacing);
    derivRow = C(2, :);
    for i = (halfWin+1):(n-halfWin)
        dy(i) = derivRow * y(i-halfWin:i+halfWin);
    end
end

%% ========================================================================
function runs = findRuns(mask)
%FINDRUNS Nx2 matrix of [startIdx, endIdx] for contiguous true runs.

    runs  = zeros(0, 2);
    n     = numel(mask);
    inRun = false;
    s = 1;
    for i = 1:n
        if mask(i) && ~inRun
            s = i; inRun = true;
        elseif ~mask(i) && inRun
            runs(end+1, :) = [s, i-1]; %#ok<AGROW>
            inRun = false;
        end
    end
    if inRun
        runs(end+1, :) = [s, n]; %#ok<AGROW>
    end
end

%% ========================================================================
function merged = mergeRuns(runs, d, maxGapNm)
%MERGERUNS Bridge consecutive runs separated by a small gap in d.

    if isempty(runs)
        merged = runs;
        return
    end
    merged = runs(1, :);
    for k = 2:size(runs, 1)
        gap = d(runs(k, 1)) - d(merged(end, 2));
        if gap <= maxGapNm
            merged(end, 2) = runs(k, 2);
        else
            merged(end+1, :) = runs(k, :); %#ok<AGROW>
        end
    end
end

%% ========================================================================
function steps = detectSteps(d, y, g, Parameters)
%DETECTSTEPS Find flat conductance plateaus ("steps") in one trace.
% Returns a struct array with fields avgG [G0], startPos [nm], length [nm].

    steps = struct('avgG', {}, 'startPos', {}, 'length', {});
    n = numel(d);
    if n < 10
        return
    end

    dspacing = median(diff(d));
    if ~(dspacing > 0)
        return
    end

    halfWin = max(3, round((Parameters.SGWindowNm / dspacing) / 2));
    dy = sgDerivative(y, dspacing, halfWin, Parameters.SGPolyOrder);

    quietMask = ~isnan(dy) & (abs(dy) < Parameters.DerivThresh) & ...
                (y > log10(Parameters.StepGmin)) & ...
                (y < log10(Parameters.StepGmax));

    runs = findRuns(quietMask);
    runs = mergeRuns(runs, d, Parameters.MergeGapNm);

    for k = 1:size(runs, 1)
        s = runs(k, 1); e = runs(k, 2);
        len = d(e) - d(s);
        if len >= Parameters.StepLenMin && len <= Parameters.StepLenMax
            steps(end+1).avgG   = mean(g(s:e)); %#ok<AGROW>
            steps(end).startPos = d(s);
            steps(end).length   = len;
        end
    end
end

%% ========================================================================
function data = ensureTraceLoaded(data, idx)
%ENSURETRACELOADED Populate data.cache{idx} and data.stepCache{idx} if
% they haven't been computed yet.

    if isempty(data.cache{idx})
        [d, y, g, ok] = readBreakingTrace(data.fileList{idx}, data.speed_const, data.Parameters);
        if ok
            data.cache{idx} = [d, y, g];
        else
            data.cache{idx} = zeros(0, 3);
        end
    end

    if isempty(data.stepCache{idx})
        dy_g = data.cache{idx};
        if isempty(dy_g)
            data.stepCache{idx} = struct('avgG', {}, 'startPos', {}, 'length', {});
        else
            data.stepCache{idx} = detectSteps(dy_g(:,1), dy_g(:,2), dy_g(:,3), data.Parameters);
        end
    end
end

%% ========================================================================
function keyPress(fig, evt)
    data = guidata(fig);

    switch evt.Key
        case {'rightarrow', 'd'}
            data.currentIndex = min(data.currentIndex + 1, data.Nfiles);
        case {'leftarrow', 'a'}
            data.currentIndex = max(data.currentIndex - 1, 1);
        case 'space'
            data.excluded(data.currentIndex) = ~data.excluded(data.currentIndex);
        case 'c'
            data.excluded(:) = false;
        case 'g'
            answer = inputdlg('Jump to trace #:', 'Go to trace', 1, {num2str(data.currentIndex)});
            if ~isempty(answer)
                n = round(str2double(answer{1}));
                if ~isnan(n)
                    data.currentIndex = min(max(n, 1), data.Nfiles);
                end
            end
        case 'e'
            guidata(fig, data);
            exportSteps(fig);
            return
        case 'q'
            close(fig);
            return
        otherwise
            return
    end

    guidata(fig, data);
    renderCurrent(fig);
end

%% ========================================================================
function renderCurrent(fig)
    data = guidata(fig);
    idx  = data.currentIndex;

    data = ensureTraceLoaded(data, idx);
    guidata(fig, data);

    dyg   = data.cache{idx};
    steps = data.stepCache{idx};
    ax    = data.ax;
    cla(ax);
    hold(ax, 'on');

    if isempty(dyg)
        text(ax, 0.5, 0.5, 'Could not read this trace', ...
             'HorizontalAlignment', 'center', 'Units', 'normalized');
    else
        d = dyg(:,1); y = dyg(:,2);
        ylo = min(y) - 0.5;
        yhi = max(y) + 0.5;

        % Highlight detected steps first, so the trace draws on top
        for k = 1:numel(steps)
            s = steps(k).startPos;
            e = s + steps(k).length;
            fill(ax, [s e e s], [ylo ylo yhi yhi], [1 0.62 0.2], ...
                 'FaceAlpha', 0.30, 'EdgeColor', 'none');
            text(ax, (s+e)/2, yhi - 0.15, sprintf('#%d', k), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.6 0.35 0]);
        end

        plot(ax, d, y, 'b-', 'LineWidth', 1);
        ylim(ax, [ylo, yhi]);
    end

    set(ax, 'Color', 'w');
    xlabel(ax, 'Displacement (nm)');
    ylabel(ax, 'log_{10}(G/G_0)');
    grid(ax, 'on');
    box(ax, 'on');
    hold(ax, 'off');

    [~, fname] = fileparts(data.fileList{idx});
    exTag = '';
    if data.excluded(idx)
        exTag = '   \color{red}[EXCLUDED FROM EXPORT]';
    end
    title(ax, sprintf('Trace %d / %d   -   %s   (%d steps)%s', ...
          idx, data.Nfiles, fname, numel(steps), exTag), 'Interpreter', 'tex');

    if ~isempty(steps)
        fprintf('Trace %d (%s): %d step(s)\n', idx, fname, numel(steps));
        for k = 1:numel(steps)
            fprintf('   #%d  avgG = %.4e G0   start = %.4f nm   length = %.4f nm\n', ...
                    k, steps(k).avgG, steps(k).startPos, steps(k).length);
        end
    end

    nExcl = sum(data.excluded);
    set(fig, 'Name', sprintf(['Step Finder Navigator  |  %d excluded  |  ' ...
                         'Arrows: browse | Space: exclude | C: clear | ' ...
                         'G: go to # | E: export | Q: quit'], nExcl));
    drawnow;
end

%% ========================================================================
function exportSteps(fig)
    data = guidata(fig);
    P    = data.Parameters;

    minK = P.ExportMinSteps;
    maxK = P.ExportMaxSteps;
    groups = cell(maxK - minK + 1, 1);   % groups{k-minK+1} = cell rows

    fprintf('\nExporting steps for %d trace(s)... (this may take a moment)\n', data.Nfiles);

    nSkippedExcluded = 0;
    nSkippedRange    = 0;
    nUnreadable      = 0;

    for idx = 1:data.Nfiles
        data = ensureTraceLoaded(data, idx);

        if data.excluded(idx)
            nSkippedExcluded = nSkippedExcluded + 1;
            continue
        end
        if isempty(data.cache{idx})
            nUnreadable = nUnreadable + 1;
            continue
        end

        steps = data.stepCache{idx};
        nSteps = numel(steps);
        if nSteps < minK || nSteps > maxK
            nSkippedRange = nSkippedRange + 1;
            continue
        end

        [~, fname] = fileparts(data.fileList{idx});
        gIdx = nSteps - minK + 1;
        for k = 1:nSteps
            groups{gIdx}(end+1, :) = {fname, steps(k).avgG, steps(k).startPos, steps(k).length}; %#ok<AGROW>
        end
    end

    guidata(fig, data);  % persist any newly-computed cache/stepCache

    outFolder = fullfile(data.folderPath, 'StepAnalysis');
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    filesWritten = {};
    for gIdx = 1:numel(groups)
        rows = groups{gIdx};
        if isempty(rows)
            continue
        end
        k = minK + gIdx - 1;
        outFile = fullfile(outFolder, sprintf('%s_%dsteps.txt', data.prefix, k));

        fid = fopen(outFile, 'w');
        fprintf(fid, 'TraceFile\tStepAvgG_G0\tStepStart_nm\tStepLength_nm\n');
        for r = 1:size(rows, 1)
            fprintf(fid, '%s\t%.6e\t%.4f\t%.4f\n', rows{r,1}, rows{r,2}, rows{r,3}, rows{r,4});
        end
        fclose(fid);

        filesWritten{end+1} = outFile; %#ok<AGROW>
        fprintf('  wrote %s  (%d step-rows from %d trace(s))\n', outFile, size(rows,1), numel(unique(rows(:,1))));
    end

    fprintf('Done. %d file(s) written to %s\n', numel(filesWritten), outFolder);
    fprintf('  skipped (excluded): %d, skipped (step count outside [%d,%d]): %d, unreadable: %d\n', ...
            nSkippedExcluded, minK, maxK, nSkippedRange, nUnreadable);

    if isempty(filesWritten)
        msgbox('No traces matched the export criteria - no files were written. Check ExportMinSteps/ExportMaxSteps and your exclusions.', 'Export complete');
    else
        msgbox(sprintf('Exported %d file(s) to:\n%s', numel(filesWritten), outFolder), 'Export complete');
    end
end