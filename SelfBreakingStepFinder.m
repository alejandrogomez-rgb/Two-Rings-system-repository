function SelfBreakingStepFinder(startFile)
%SELFBREAKINGSTEPFINDER Interactive browser + detector for conductance
% steps (plateaus) in self-breaking MCBJ traces.
%
% This is the self-breaking counterpart of StepFinderNavigator.m (which
% was built for fast-breaking traces on a piezo-displacement axis), using
% a ROLLING STANDARD DEVIATION of y = log10(G/G0) instead of a
% Savitzky-Golay derivative, on a TIME [s] axis instead of displacement
% [nm] (self-breaking has no piezo ramp, so there is no d = f(t)
% conversion).
%
% Unlike the original flat-plateau finder, the rolling STD here is used
% to find TRANSITIONS (jumps between conductance states), not steps
% directly. A "step" is then simply any quiet state bounded by a
% transition -- whether its average conductance is higher or lower than
% that of the quiet state immediately preceding it in time no longer
% matters.
%
% ---- How a step is detected -------------------------------------------
%   1. t [s] and y = log10(G/G0) are read directly from the file (no unit
%      conversion needed -- self-breaking traces are already time-based).
%   2. A rolling standard deviation of y is computed over a sliding
%      window of width Parameters.StepWindowSec.
%   3. Points where rollingSTD < Parameters.TransitionStdMin are flagged
%      as "quiet" (part of a stable state, not a transition). This is
%      STD-only: Parameters.StepGmin/StepGmax are deliberately NOT
%      applied here, so that G merely wobbling across that window's edge
%      (while the rolling STD stays low) is never mistaken for a real
%      transition.
%   4. Consecutive quiet points are grouped into runs; short interruptions
%      (noise blips that briefly push the STD over threshold) are
%      bridged if the gap is below Parameters.MergeGapSec -- so only
%      gaps wider than MergeGapSec survive as genuine transitions.
%   5. Every run (candidate state) is checked in time order and kept as
%      a genuine "step" only if ALL of the following hold:
%        - its average conductance falls within (Parameters.StepGmin,
%          Parameters.StepGmax);
%        - it is relatively stable/flat: |y_end - y_start| across the
%          state is below Parameters.StepMaxDriftDec decades;
%        - its duration falls within [Parameters.StepLenMin,
%          Parameters.StepLenMax] seconds ("m" and "n");
%        - it is at a genuinely different conductance level than the
%          PREVIOUS ACCEPTED STEP (not just whatever run happens to
%          precede it in time -- a run rejected by any of the above
%          checks is skipped over and never used for this comparison):
%          |log10(avgG) - log10(avgG of the previous accepted step)|
%          is at least Parameters.StepMinSeparationDec decades. The
%          very first candidate has no previous step to compare
%          against, so it is exempt from this check.
%
% USAGE:
%   SelfBreakingStepFinder()            % prompts you to pick a .dat file
%   SelfBreakingStepFinder(fullpath)    % start directly from a known file
%
% Any file matching the same "<prefix>_<number>.dat" naming as the file
% you pick/pass will be included, sorted by that trailing number (this
% matches e.g. "Self_Breaking_260712_0823_2.dat" -> prefix
% "Self_Breaking_260712_0823", number 2).
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
% ---- EDIT THIS IF YOUR SETUP CHANGES ---------------------------
Parameters.minlvl = 1e-8;   % G0; conductance floor before log10 (no baseline offset)
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---- STEP-DETECTION PARAMETERS: tune these while browsing -----
Parameters.StepWindowSec   = 0.75;   % s, total width of the rolling-STD window
Parameters.TransitionStdMin = 0.2;  % decades of G0; rolling STD AT or ABOVE this marks a transition. Below it = a quiet candidate state.
Parameters.MergeGapSec     = 0.25;   % s; bridges brief noise-driven interruptions (gaps shorter than this are NOT treated as real transitions)
Parameters.StepGmin        = 1e-7;   % G0; states are only ever considered above this
Parameters.StepGmax        = 1e-1;     % G0; states are only ever considered below this
Parameters.StepLenMin      = 0.5;    % s; minimum accepted step length
Parameters.StepLenMax      = 30;     % s; maximum accepted step length
Parameters.StepMaxDriftDec = 1;    % decades; a step is rejected if |y_end - y_start| across it is AT or ABOVE this (keeps only relatively stable/flat states)
Parameters.StepMinSeparationDec = 0.3; % decades; a step is rejected unless its avgG differs from the avgG of the PREVIOUS ACCEPTED STEP (not just the preceding raw run) by AT LEAST this much (keeps consecutive steps from being the same conductance value)
%   -> Too few / fragmented states found?  Increase TransitionStdMin or MergeGapSec.
%   -> States blurring together / noise picked up as a state? Decrease them,
%      or shrink StepWindowSec so edges are less smeared.
%   -> Use the bottom panel (rolling STD vs. the red TransitionStdMin line) to
%      pick a threshold that separates stable states from real transitions.
%   -> Steps drifting too much across their own length? Decrease StepMaxDriftDec.
%   -> Getting steps that are basically the same conductance as the one
%      right before them? Increase StepMinSeparationDec.
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---- EXPORT PARAMETERS -----------------------------------------
Parameters.ExportMinSteps = 3;    % smallest step-count group to export
Parameters.ExportMaxSteps = 10;    % largest step-count group to export
% Only step counts in [ExportMinSteps, ExportMaxSteps] that have at
% least one matching (non-excluded) trace will produce a .txt file.
% ---------------------------------------------------------------

    % ---- 1) Pick a representative file, if not given ----
    if nargin < 1 || isempty(startFile)
        [fname, fpath] = uigetfile('*.dat', 'Select any .dat file from the self-breaking series');
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
    fprintf('Found %d self-breaking trace file(s) matching "%s_*.dat"\n', Nfiles, prefix);

    % ---- 3) Build the navigator figure (two linked panels) ----
    fig = figure('Color', 'w', 'Name', 'Self-Breaking Step Finder', ...
                 'NumberTitle', 'off', 'KeyPressFcn', @keyPress);
    ax1 = subplot(2, 1, 1, 'Parent', fig);
    ax2 = subplot(2, 1, 2, 'Parent', fig);
    set(ax1, 'Color', 'w');
    set(ax2, 'Color', 'w');
    linkaxes([ax1, ax2], 'x');

    data.fileList     = fileList;
    data.Nfiles        = Nfiles;
    data.currentIndex  = 1;
    data.excluded      = false(1, Nfiles);
    data.cache         = cell(Nfiles, 1);   % each entry: [t, y, g]
    data.stepCache     = cell(Nfiles, 1);   % each entry: steps struct array
    data.rstdCache     = cell(Nfiles, 1);   % each entry: rolling-STD vector
    data.Parameters    = Parameters;
    data.ax1           = ax1;
    data.ax2           = ax2;
    data.folderPath    = folderPath;
    data.prefix        = prefix;
    guidata(fig, data);

    fprintf(['\nControls: Right/d = next, Left/a = previous, Space = exclude/include, ' ...
             'c = clear exclusions, g = go to #, e = export grouped .txt, q = quit\n\n']);

    renderCurrent(fig);
end

%% ========================================================================
function [t, y, g, ok] = readSelfBreakingTrace(fullFilename, Parameters)
%READSELFBREAKINGTRACE Read the self-breaking segment (right after the
% first "@" section marker, up to the next non-numeric line -- typically
% either a second "@" marker starting a trailing verification pulse, or
% end-of-file) of one .dat file. Returns time t [s], log10(G/G0) y, and
% the (floor-clipped) conductance g [G0] used to compute y.

    t = []; y = []; g = []; ok = false;

    fid = fopen(fullFilename, 'r');
    fc  = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    fc  = fc{:};

    atLines = find(strncmp('@', fc, 1));
    if isempty(atLines)
        warning('%s: no "@" section marker found, skipping.', fullFilename);
        return
    end
    if numel(fc) <= atLines(1)
        warning('%s: no data after "@" marker, skipping.', fullFilename);
        return
    end

    % Detect number of columns from the first data line
    tmp  = regexp(fc(atLines(1)+1), '\s*([\S]*)\s*', 'tokens');
    cols = numel(tmp{1});

    fid = fopen(fullFilename, 'r');
    if cols <= 2
        raw = textscan(fid, '%f %f', 'HeaderLines', atLines(1));
    else
        raw = textscan(fid, '%f %f %f', 'HeaderLines', atLines(1));
    end
    fclose(fid);

    t = raw{1};
    g = raw{2};

    if isempty(t)
        warning('%s: empty data section, skipping.', fullFilename);
        return
    end

    t = t - t(1);   % zero-reference time (usually already starts at 0)

    g(g < Parameters.minlvl) = Parameters.minlvl;
    y  = log10(g);
    ok = true;
end

%% ========================================================================
function rs = rollingStd(y, halfWin)
%ROLLINGSTD Dependency-free rolling standard deviation (no toolbox / no
% movstd needed -- keeps this Octave-compatible). Edge points (within
% halfWin of either end) are returned as NaN.

    n  = numel(y);
    rs = nan(n, 1);
    if n < (2*halfWin + 1)
        return
    end
    for i = (halfWin+1):(n-halfWin)
        rs(i) = std(y(i-halfWin:i+halfWin));
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
function merged = mergeRuns(runs, t, maxGapSec)
%MERGERUNS Bridge consecutive runs separated by a small gap in time.

    if isempty(runs)
        merged = runs;
        return
    end
    merged = runs(1, :);
    for k = 2:size(runs, 1)
        gap = t(runs(k, 1)) - t(merged(end, 2));
        if gap <= maxGapSec
            merged(end, 2) = runs(k, 2);
        else
            merged(end+1, :) = runs(k, :); %#ok<AGROW>
        end
    end
end

%% ========================================================================
function [steps, rstd] = detectSteps(t, y, g, Parameters)
%DETECTSTEPS Find conductance steps in one self-breaking trace by first
% finding TRANSITIONS (high rolling-STD regions) and then treating every
% quiet state flanked by a transition as a candidate step. A candidate is
% kept only if it also falls within the G/duration windows, is
% relatively flat (small drift between its start and end), and sits at a
% conductance level clearly separated from the PREVIOUS ACCEPTED STEP
% (rejected runs are skipped over and don't count). Returns a struct
% array with fields avgG [G0],
% startTime [s], length [s], stdY (STD of y within the step, for
% reference), plus the full rolling-STD vector (for plotting/tuning).

    steps = struct('avgG', {}, 'startTime', {}, 'length', {}, 'stdY', {});
    n = numel(t);
    rstd = nan(n, 1);
    if n < 10
        return
    end

    dt = median(diff(t));
    if ~(dt > 0)
        return
    end

    halfWin = max(2, round((Parameters.StepWindowSec / dt) / 2));
    rstd = rollingStd(y, halfWin);

    % Quiet points = candidate states (NOT yet "steps"); everything else
    % (rolling STD AT or ABOVE TransitionStdMin) is a candidate
    % transition. NOTE: StepGmin/StepGmax are intentionally NOT applied
    % here -- see point 3 above. They're applied once per COMPLETE
    % candidate state below instead.
    quietMask = ~isnan(rstd) & (rstd < Parameters.TransitionStdMin);

    runs = findRuns(quietMask);
    runs = mergeRuns(runs, t, Parameters.MergeGapSec);
    % Each row of runs is now one candidate state; the gap between
    % consecutive rows is one surviving transition (gaps shorter than
    % MergeGapSec were already bridged away above).

    nRuns = size(runs, 1);
    if nRuns < 2
        return   % need at least one transition to validate anything
    end

    avgGs = zeros(nRuns, 1);
    for k = 1:nRuns
        avgGs(k) = mean(g(runs(k,1):runs(k,2)));
    end

    % Each run is checked against the running record of the last
    % ACCEPTED step (not merely whatever run happens to precede it in
    % time) -- so a run that got rejected for being out of range, too
    % short, or too drifty is skipped over and does NOT count as "the
    % previous step" for separation purposes. lastStepG stays empty
    % until the first step is accepted, so that first candidate is never
    % rejected on separation grounds.
    lastStepG = [];
    for k = 1:nRuns
        if avgGs(k) <= Parameters.StepGmin || avgGs(k) >= Parameters.StepGmax
            continue
        end
        s = runs(k, 1); e = runs(k, 2);
        if abs(y(e) - y(s)) >= Parameters.StepMaxDriftDec
            continue
        end
        len = t(e) - t(s);
        if ~(len >= Parameters.StepLenMin && len <= Parameters.StepLenMax)
            continue
        end
        if ~isempty(lastStepG)
            sepDec = abs(log10(avgGs(k)) - log10(lastStepG));
            if sepDec < Parameters.StepMinSeparationDec
                continue
            end
        end

        steps(end+1).avgG      = avgGs(k); %#ok<AGROW>
        steps(end).startTime   = t(s);
        steps(end).length      = len;
        steps(end).stdY        = std(y(s:e));
        lastStepG = avgGs(k);
    end
end

%% ========================================================================
function data = ensureTraceLoaded(data, idx)
%ENSURETRACELOADED Populate data.cache{idx}, data.stepCache{idx} and
% data.rstdCache{idx} if they haven't been computed yet.

    if isempty(data.cache{idx})
        [t, y, g, ok] = readSelfBreakingTrace(data.fileList{idx}, data.Parameters);
        if ok
            data.cache{idx} = [t, y, g];
        else
            data.cache{idx} = zeros(0, 3);
        end
    end

    if isempty(data.stepCache{idx})
        tyg = data.cache{idx};
        if isempty(tyg)
            data.stepCache{idx} = struct('avgG', {}, 'startTime', {}, 'length', {}, 'stdY', {});
            data.rstdCache{idx} = [];
        else
            [steps, rstd] = detectSteps(tyg(:,1), tyg(:,2), tyg(:,3), data.Parameters);
            data.stepCache{idx} = steps;
            data.rstdCache{idx} = rstd;
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
    P    = data.Parameters;

    data = ensureTraceLoaded(data, idx);
    guidata(fig, data);

    tyg   = data.cache{idx};
    steps = data.stepCache{idx};
    rstd  = data.rstdCache{idx};
    ax1   = data.ax1;
    ax2   = data.ax2;

    cla(ax1); hold(ax1, 'on');
    cla(ax2); hold(ax2, 'on');

    if isempty(tyg)
        text(ax1, 0.5, 0.5, 'Could not read this trace', ...
             'HorizontalAlignment', 'center', 'Units', 'normalized');
    else
        t = tyg(:,1); y = tyg(:,2);
        xlim(ax1, [t(1) t(end)]);
        ylo = min(y) - 0.5;
        yhi = max(y) + 0.5;

        rlo = 0;
        if isempty(rstd) || all(isnan(rstd))
            rhi = max(P.TransitionStdMin * 2, 0.1);
        else
            rhi = max(P.TransitionStdMin * 1.5, max(rstd(~isnan(rstd))) * 1.05);
        end

        % Highlight detected steps first, so the traces draw on top
        for k = 1:numel(steps)
            s = steps(k).startTime;
            e = s + steps(k).length;
            fill(ax1, [s e e s], [ylo ylo yhi yhi], [1 0.62 0.2], ...
                 'FaceAlpha', 0.30, 'EdgeColor', 'none');
            text(ax1, (s+e)/2, yhi - 0.15, sprintf('#%d', k), ...
                 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.6 0.35 0]);
            fill(ax2, [s e e s], [rlo rlo rhi rhi], [1 0.62 0.2], ...
                 'FaceAlpha', 0.30, 'EdgeColor', 'none');
        end

        % Reference lines for the conductance window on the top panel
        plot(ax1, [t(1) t(end)], log10(P.StepGmin)*[1 1], 'k--', 'LineWidth', 0.75);
        plot(ax1, [t(1) t(end)], log10(P.StepGmax)*[1 1], 'k--', 'LineWidth', 0.75);

        plot(ax1, t, y, 'b-', 'LineWidth', 1);
        ylim(ax1, [ylo, yhi]);

        % Rolling STD + threshold line on the bottom panel
        plot(ax2, t, rstd, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1);
        plot(ax2, [t(1) t(end)], P.TransitionStdMin*[1 1], 'r--', 'LineWidth', 1);
        ylim(ax2, [rlo, rhi]);
    end

    set(ax1, 'Color', 'w');
    ylabel(ax1, 'log_{10}(G/G_0)');
    grid(ax1, 'on'); box(ax1, 'on'); hold(ax1, 'off');

    set(ax2, 'Color', 'w');
    xlabel(ax2, 'Time (s)');
    ylabel(ax2, 'Rolling STD (dec)');
    grid(ax2, 'on'); box(ax2, 'on'); hold(ax2, 'off');

    [~, fname] = fileparts(data.fileList{idx});
    exTag = '';
    if data.excluded(idx)
        exTag = '   \color{red}[EXCLUDED FROM EXPORT]';
    end
    title(ax1, sprintf('Trace %d / %d   -   %s   (%d steps)%s', ...
          idx, data.Nfiles, fname, numel(steps), exTag), 'Interpreter', 'tex');

    if ~isempty(steps)
        fprintf('Trace %d (%s): %d step(s)\n', idx, fname, numel(steps));
        for k = 1:numel(steps)
            fprintf('   #%d  avgG = %.4e G0   start = %.4f s   length = %.4f s   stdY = %.4f\n', ...
                    k, steps(k).avgG, steps(k).startTime, steps(k).length, steps(k).stdY);
        end
    end

    nExcl = sum(data.excluded);
    set(fig, 'Name', sprintf(['Self-Breaking Step Finder  |  %d excluded  |  ' ...
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
            groups{gIdx}(end+1, :) = {fname, steps(k).avgG, steps(k).startTime, steps(k).length}; %#ok<AGROW>
        end
    end

    guidata(fig, data);  % persist any newly-computed cache/stepCache/rstdCache

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
        fprintf(fid, 'TraceFile\tStepAvgG_G0\tStepStartTime_s\tStepLength_s\n');
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
