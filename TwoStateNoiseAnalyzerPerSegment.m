function TwoStateNoiseAnalyzerPerSegment(startFile)
%TWOSTATENOISEANALYZERPERSEGMENT Interactive browser + two-state
% (switching vs. stable) classifier for self-breaking MCBJ traces, with
% PER-SEGMENT export.
%
% This is a variant of TwoStateNoiseAnalyzer.m. Trace loading, trimming,
% windowed-std discrimination, segmentation, gap-merging and the
% trace-level "both states present" gate all work IDENTICALLY to that
% script. The difference is entirely in what gets exported:
%   - TwoStateNoiseAnalyzer.m combines all segments sharing a label within
%     a trace into ONE aggregated row per state per trace (duration-
%     weighted-averaged G and std).
%   - TwoStateNoiseAnalyzerPerSegment.m instead writes ONE ROW PER
%     INDIVIDUAL SEGMENT - no combining/averaging across segments of the
%     same state. A trace with, say, 3 switching segments and 2 stable
%     segments contributes 5 rows, not 2.
% The 2D histogram is built to match: it now has one count PER SEGMENT
% (instead of one count per state per trace as in the original script).
% An optional duration-weighting mode is also available (see
% Parameters.WeightBinsBySegmentDuration below): when enabled, each
% segment contributes its own duration to its histogram bin instead of a
% flat count of 1, so longer segments carry proportionally more weight.
%
% The trace-level "both states present" gate is UNCHANGED from
% TwoStateNoiseAnalyzer.m: a trace is only exported if the COMBINED
% (aggregated) length of each state exceeds Parameters.MinStateLengthS,
% exactly as before. A separate, output-only parameter,
% Parameters.MinSegmentDurationOutputS, additionally prunes which
% INDIVIDUAL segments get written out/histogrammed once a trace has
% already passed that gate - it plays NO role in deciding whether a trace
% qualifies in the first place.
%
% ---- How classification works ------------------------------------------
%   1. d [s] (zero-referenced time, used as the x-axis), t [s] (raw time)
%      and y = log10(G/G0) are computed in readBreakingTrace. Self-breaking
%      traces are recorded at fixed bias, i.e. G vs time, so d is
%      just t shifted to zero at the onset of breaking - there is no
%      piezo-driven displacement.
%   2. The trace is trimmed to the window between the first point where
%      G drops to/below Parameters.TrimStartG and the first (subsequent)
%      point where G drops to/below Parameters.TrimEndG. That start point
%      is then pushed forward by a further fixed Parameters.TrimStartTimeS
%      seconds, and that end point is pulled back by a further fixed
%      Parameters.TrimEndTimeS seconds, to cut out high-std regions
%      sitting right after/before those G-based trim points.
%   3. A windowed standard-deviation curve r(i) = std(y) [decades] is
%      computed over the trimmed region, using a centered window of width
%      WindowSizeS.
%   4. Segments are formed by thresholding r(i) directly against
%      Parameters.StdThresh and grouping contiguous runs.
%   5. Any segment that is sandwiched between two same-label neighbours
%      (i.e. a short interruption of the opposite state) is merged into
%      those neighbours if its duration is below Parameters.MaxMergeGapS.
%   6. For the TRACE GATE ONLY: all segments sharing the same label
%      ("switching"/"stable") within one trace are combined into a single
%      aggregated state (durations summed, G duration-weighted-averaged).
%      A state only counts as "present" in the trace if its combined
%      length exceeds Parameters.MinStateLengthS (this is deliberately
%      applied to the COMBINED length, not each individual fragment -
%      simple thresholding can still chop one real switching region into
%      several short fragments if the signal briefly dips near the
%      threshold, and requiring each fragment alone to clear the minimum
%      would wrongly discard genuine, just-fragmented, switching regions).
%      Only traces where BOTH states are "present" this way are exported.
%      This aggregation is used SOLELY to decide the gate - the actual
%      exported rows are the individual segments themselves (see above),
%      each optionally filtered by Parameters.MinSegmentDurationOutputS.
%
% USAGE:
%   TwoStateNoiseAnalyzerPerSegment()            % prompts you to pick a .dat file
%   TwoStateNoiseAnalyzerPerSegment(fullpath)    % start directly from a known file
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
%   e                 -> export .txt files for ALL traces
%   q                 -> close the navigator
%
% ---------------------------------------------------------------
% ---- EDIT THESE IF YOUR SETUP CHANGES (same as loadRawData.m) ----
Parameters.ZeroRefG       = 0.3;        % G/G0 threshold used to zero-reference the time axis
                                         % (t=0 set at the first point where G drops below this)
Parameters.minlvl         = 1e-8;       % conductance floor before log10 (no baseline offset)
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---- TWO-STATE DETECTION PARAMETERS: tune these while browsing ------
Parameters.WindowSizeS     = 0.15;    % seconds, width of the
                                      % sliding discriminator window
Parameters.StdThresh       = 0.1;   % decades of G0; window std above this = "switching"
Parameters.TrimStartG      = 1e-1;   % G0, analysis window starts at the first point where
                                      % G drops to/below this (everything before is trimmed)
Parameters.TrimEndG        = 1e-7;   % G0, analysis window ends at the first point (at or
                                      % after TrimStartG) where G drops to/below this
                                      % (everything after is trimmed)
Parameters.TrimStartTimeS  = 1;      % seconds, ADDITIONAL trim
                                      % applied AFTER the G-based start trim above: pushes
                                      % dStart forward by this much more time, to cut out the
                                      % high-std region that tends to sit right after the
                                      % G-based start trim point. Set to 0 to disable.
Parameters.TrimEndTimeS    = 1;    % seconds, ADDITIONAL trim
                                      % applied AFTER the G-based end trim above: shortens
                                      % dEnd by this much more time, to cut out the high-std
                                      % region that tends to sit right before the G-based
                                      % trim point. Set to 0 to disable.
Parameters.MaxMergeGapS    = 0.15;    % seconds; a segment that is
                                      % sandwiched between two same-label neighbours (i.e. a
                                      % short interruption of the opposite state) is merged
                                      % into those neighbours - and the merged run takes on
                                      % their label - if its duration is below this value.
                                      % Set to 0 to disable merging.
Parameters.MinStateLengthS = 2.5;       % seconds; a state's COMBINED
                                      % length (summed over all its segments in the trace) must
                                      % exceed this to count as "present" - filters out spurious
                                      % single-window blips. USED ONLY FOR THE TRACE-LEVEL "both
                                      % states present" GATE - see Parameters.MinSegmentDurationOutputS
                                      % below for the (separate) per-segment OUTPUT filter.
%   -> Missing an obvious switching region? Lower StdThresh.
%   -> Picking up switching everywhere / too sensitive? Raise StdThresh, or
%      raise MinStateLengthS to ignore short spurious detections.
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---- HISTOGRAM PARAMETERS (for the segment_histogram2d_per_segment.txt export) ----
Parameters.HistNBinsG    = 30;   % number of log-spaced bins along mean-G
Parameters.HistNBinsStat = 30;   % number of linear bins along std
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% ---- PER-SEGMENT EXPORT PARAMETERS ----
Parameters.MinSegmentDurationOutputS = 0.25;   % seconds; INDIVIDUAL segments
                                      % shorter than this (by durationS) are left out of the
                                      % per-segment export (summary/scatter/histogram) ONLY.
                                      % This does NOT feed back into the trace-level "both
                                      % states present" gate above (that gate always uses the
                                      % full, unfiltered aggregate vs. Parameters.MinStateLengthS,
                                      % exactly as in TwoStateNoiseAnalyzer.m). Set to 0 to
                                      % disable (export every segment regardless of length).
Parameters.WeightBinsBySegmentDuration = true;   % boolean; if true, each
                                      % segment contributes its own durationS (in seconds) to
                                      % its 2D-histogram bin instead of a flat count of 1, so
                                      % longer segments have a proportionally bigger effect on
                                      % the histogram. If false (default), every exported
                                      % segment contributes an unweighted count of 1, matching
                                      % the counting convention of TwoStateNoiseAnalyzer.m
                                      % (just at per-segment instead of per-state resolution).
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

    % ---- 3) Build the navigator figure ----
    fig = figure('Color', 'w', 'Name', 'Two-State Noise Analyzer (Per-Segment)', ...
                 'NumberTitle', 'off', 'KeyPressFcn', @keyPress);
    ax = axes(fig);
    set(ax, 'Color', 'w');

    data.fileList     = fileList;
    data.Nfiles        = Nfiles;
    data.currentIndex  = 1;
    data.excluded      = false(1, Nfiles);
    data.cache         = cell(Nfiles, 1);   % each entry: struct with d,t,y,g
    data.segCache      = cell(Nfiles, 1);   % each entry: segs struct array
    data.aggCache      = cell(Nfiles, 1);   % each entry: agg struct (switching/stable)
    data.Parameters    = Parameters;
    data.ax            = ax;
    data.folderPath    = folderPath;
    data.prefix        = prefix;
    guidata(fig, data);

    fprintf(['\nControls: Right/d = next, Left/a = previous, Space = exclude/include, ' ...
             'c = clear exclusions, g = go to #, e = export, q = quit\n\n']);

    renderCurrent(fig);
end

%% ========================================================================
function [d, t, y, g, ok] = readBreakingTrace(fullFilename, Parameters)
%READBREAKINGTRACE Read the breaking segment (between the 1st and 2nd "@")
% of one .dat file and return zero-referenced time d [s] (used as the
% x-axis - self-breaking traces are recorded at fixed bias, i.e. G vs
% time, so there is no piezo displacement to convert to here), raw time
% t [s] (NOT zero-referenced), log10(G/G0) y, and the (floor-clipped)
% conductance g [G0] used to compute y.

    d = []; t = []; y = []; g = []; ok = false;

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

    % Zero-reference time at first point where g < ZeroRefG, so the x-axis
    % (d) reads 0 s at the onset of breaking, consistently across traces.
    d = t;   % seconds
    try
        idx0   = find(g < Parameters.ZeroRefG);
        idx0   = idx0(1);
        offset = d(idx0);
        d      = d - offset;
    catch
        % trace never drops below ZeroRefG -> leave d unshifted
    end

    g(g < Parameters.minlvl) = Parameters.minlvl;
    y  = log10(g);
    ok = true;
end

%% ========================================================================
function r = rollingStd(y, halfWin)
%ROLLINGSTD Windowed standard deviation of y, centered window of
% (2*halfWin+1) samples. Edge points (within halfWin of either end) are
% NaN. Vectorized via cumulative sums for speed (O(n) total).

    n = numel(y);
    r = nan(n, 1);
    if n < (2*halfWin + 1)
        return
    end
    cs  = [0; cumsum(y)];
    cs2 = [0; cumsum(y.^2)];
    for i = (halfWin+1):(n-halfWin)
        a = i - halfWin; b = i + halfWin; N = b - a + 1;
        s1 = cs(b+1) - cs(a);
        s2 = cs2(b+1) - cs2(a);
        m  = s1 / N;
        v  = max(s2/N - m^2, 0);
        r(i) = sqrt(v * N / (N-1));
    end
end

%% ========================================================================
function segs = groupRunsByLabel(labelVec)
%GROUPRUNSBYLABEL Group a vector of {1, 0, NaN} into contiguous runs of
% the same non-NaN value. NaN samples break a run but do not themselves
% start one. Returns Nx3 [startIdx, endIdx, label].

    segs = zeros(0, 3);
    n = numel(labelVec);
    curLabel = NaN;
    s = 0;
    for i = 1:n
        v = labelVec(i);
        if isnan(v)
            if ~isnan(curLabel)
                segs(end+1, :) = [s, i-1, curLabel]; %#ok<AGROW>
                curLabel = NaN;
            end
            continue
        end
        if isnan(curLabel)
            curLabel = v; s = i;
        elseif v ~= curLabel
            segs(end+1, :) = [s, i-1, curLabel]; %#ok<AGROW>
            curLabel = v; s = i;
        end
    end
    if ~isnan(curLabel)
        segs(end+1, :) = [s, n, curLabel]; %#ok<AGROW>
    end
end

%% ========================================================================
function runList = mergeShortGaps(runList, d, maxGapS)
%MERGESHORTGAPS Fold short "gap" runs back into the surrounding state.
% runList is the Nx3 [startIdx, endIdx, label] array produced by
% groupRunsByLabel (assumed to already be in index order, with each run's
% endIdx immediately preceding the next run's startIdx). Any run whose
% label differs from BOTH its immediate neighbours (i.e. it's a short
% interruption of the opposite state, sandwiched between two same-label
% runs) is absorbed into those neighbours - all three become one run
% carrying the neighbours' label - provided the gap run's duration
% (measured on d, the zero-referenced time axis) is below maxGapS.
% Repeats until no more qualifying gaps remain, since merging can expose
% new same-label neighbour pairs.

    if maxGapS <= 0 || isempty(runList) || size(runList, 1) < 3
        return
    end

    changed = true;
    while changed
        changed = false;
        k = 2;
        while k <= size(runList, 1) - 1
            prevLabel = runList(k-1, 3);
            curLabel  = runList(k,   3);
            nextLabel = runList(k+1, 3);
            if curLabel ~= prevLabel && prevLabel == nextLabel
                gapDurS = d(runList(k, 2)) - d(runList(k, 1));
                if gapDurS < maxGapS
                    runList(k-1, 2) = runList(k+1, 2);   % extend prev run to swallow gap + next
                    runList(k+1, :) = [];
                    runList(k,   :) = [];
                    changed = true;
                    break   % restart the scan since indices shifted
                end
            end
            k = k + 1;
        end
    end
end

%% ========================================================================
function [dStart, dEnd] = conductanceTrimBounds(d, g, Parameters)
%CONDUCTANCETRIMBOUNDS Find the [dStart, dEnd] analysis window (in the
% same units as d) bounded by conductance thresholds instead of fixed
% time amounts: dStart is where G first drops to/below Parameters.TrimStartG,
% dEnd is the first point at/after that where G drops to/below
% Parameters.TrimEndG. Falls back to the full trace extent on either end
% if the corresponding threshold is never crossed. dStart is then pushed
% forward by Parameters.TrimStartTimeS and dEnd pulled back further by
% Parameters.TrimEndTimeS (fixed time cuts applied AFTER the G-based
% trim, to exclude high-std regions that tend to sit right after/before
% the G-based trim points). Both additional trims are computed against
% the original G-based bounds and then clamped so the window never
% inverts (collapsing to a single point if the two trims would overlap).

    n = numel(d);
    idxStart = find(g <= Parameters.TrimStartG, 1, 'first');
    if isempty(idxStart)
        idxStart = 1;
    end
    idxEndRel = find(g(idxStart:end) <= Parameters.TrimEndG, 1, 'first');
    if isempty(idxEndRel)
        idxEnd = n;
    else
        idxEnd = idxStart + idxEndRel - 1;
    end
    dStart0 = d(idxStart);
    dEnd0   = d(idxEnd);

    % Additional time-based trims, applied AFTER the G-based trim above.
    dStart = min(dEnd0,   dStart0 + Parameters.TrimStartTimeS);
    dEnd   = max(dStart0, dEnd0   - Parameters.TrimEndTimeS);
    if dStart > dEnd
        dMid   = (dStart + dEnd) / 2;
        dStart = dMid;
        dEnd   = dMid;
    end
end

%% ========================================================================
function segs = segmentTrace(d, t, y, g, Parameters)
%SEGMENTTRACE Trim the trace, compute the windowed-std discriminator
% curve, cut it into segments by per-window thresholding, optionally
% merge short opposite-label interruptions back into their surrounding
% state (Parameters.MaxMergeGapS), and label each resulting segment
% "switching" or "stable". Returns a struct array with fields: label,
% startIdx, endIdx, statVal, avgG, maxG, minG, durationS, lengthS.
% Operates on sample indices into the FULL (untrimmed) d/t/y/g arrays.

    segs = struct('label', {}, 'startIdx', {}, 'endIdx', {}, ...
                   'statVal', {}, 'avgG', {}, 'maxG', {}, 'minG', {}, ...
                   'durationS', {}, 'lengthS', {});

    n = numel(d);
    if n < 10
        return
    end

    dspacing = median(diff(d));
    if ~(dspacing > 0)
        return
    end
    halfWin = max(3, round((Parameters.WindowSizeS / dspacing) / 2));

    r = rollingStd(y, halfWin);
    Thresh = Parameters.StdThresh;

    % Restrict to the conductance-bounded analysis window
    [dStart, dEnd] = conductanceTrimBounds(d, g, Parameters);
    inTrim = (d >= dStart) & (d <= dEnd);
    r(~inTrim) = NaN;

    if ~any(~isnan(r))
        return
    end

    % ---- Find raw candidate segments [startIdx, endIdx] via thresholding ----
    labelVec = nan(n, 1);
    labelVec(~isnan(r)) = r(~isnan(r)) > Thresh;
    runList = groupRunsByLabel(labelVec);

    % ---- Merge short opposite-label interruptions back into their
    %      surrounding state ----
    runList = mergeShortGaps(runList, d, Parameters.MaxMergeGapS);

    % ---- Convert to output segs struct ----
    for k = 1:size(runList, 1)
        s = runList(k,1); e = runList(k,2); lab = runList(k,3);
        segs(end+1).label     = ternary(lab==1, 'switching', 'stable'); %#ok<AGROW>
        segs(end).startIdx    = s;
        segs(end).endIdx      = e;
        segs(end).statVal     = mean(r(s:e), 'omitnan');
        segs(end).avgG        = mean(g(s:e));
        segs(end).maxG        = max(g(s:e));
        segs(end).minG        = min(g(s:e));
        segs(end).durationS   = t(e) - t(s);
        segs(end).lengthS     = d(e) - d(s);
    end
end

%% ========================================================================
function out = ternary(cond, a, b)
%TERNARY Small helper: cond ? a : b (MATLAB has no ?: operator).
    if cond
        out = a;
    else
        out = b;
    end
end

%% ========================================================================
function edges = safeBinEdges(vals, nBins)
%SAFEBINEDGES Build nBins+1 linearly-spaced bin edges spanning [min(vals),
% max(vals)], guaranteed to be strictly monotonically increasing (as
% required by discretize) regardless of the data's magnitude or of ties
% among values.
%
% A plain "linspace(min(vals), max(vals) + eps, nBins+1)" (as used
% upstream in some scripts) can silently fail this: MATLAB's bare `eps`
% is the spacing between representable doubles AT 1 (~2.22e-16). For
% values whose magnitude is >> 1 (or even just not close to 1), the local
% spacing between representable doubles (eps(x)) is much larger, so
% "max(vals) + eps" can round right back down to max(vals) - collapsing
% the range to zero width and producing non-increasing edges. This is
% especially likely per-segment, where many individual segments can share
% an identical value (e.g. several "stable" segments all clipped to the
% same conductance floor, or several flat segments all with std exactly
% 0), rather than the smoothed-out, duration-weighted-averaged values
% aggregateStates used to produce per-trace.

    vals = vals(isfinite(vals));
    if isempty(vals)
        vals = 0;   % degenerate fallback: still produce a usable (if
                     % meaningless) edge vector rather than erroring
    end
    lo = min(vals);
    hi = max(vals);
    span = hi - lo;

    if span > 0
        % Pad the upper edge well past ordinary rounding error at this
        % magnitude and at this bin count, so every interior linspace
        % point remains distinct.
        pad = max(span * 1e-6, eps(hi) * (nBins + 10));
    else
        % All values identical (or numerically indistinguishable) -
        % fabricate a small span around that single value, scaled to its
        % magnitude, so binning still produces a valid (single-bin-ish)
        % result instead of erroring.
        pad = max(abs(hi), 1) * 1e-6;
    end
    hi = hi + pad;

    edges = linspace(lo, hi, nBins + 1);

    % Belt-and-suspenders: force strict monotonicity even in the
    % pathological case where floating-point rounding still produced a
    % tie somewhere in the interior of the linspace output.
    for i = 2:numel(edges)
        if edges(i) <= edges(i-1)
            edges(i) = edges(i-1) + eps(edges(i-1)) * 2;
        end
    end
end

%% ========================================================================
function agg = aggregateStates(segs, Parameters)
%AGGREGATESTATES Combine all segments of the same label into one
% aggregated state (durations/lengths summed, G and statVal
% duration-weighted-averaged). maxG/minG are NOT taken across all of the
% state's segments - they come only from the single LONGEST segment
% (by lengthS) of that state, so a brief spike/dip in a short fragment
% doesn't set the reported extremum for the whole (possibly fragmented)
% state. A state counts as "present" only if its COMBINED length exceeds
% Parameters.MinStateLengthS.
%
% NOTE: this aggregation is used ONLY to decide the trace-level "both
% states present" gate (see exportPerSegment). It plays no role in the
% per-segment rows that actually get exported - those come straight from
% the individual segments in `segs`.

    agg.switching = emptyStateAgg();
    agg.stable    = emptyStateAgg();

    for k = 1:numel(segs)
        s = segs(k);
        agg.(s.label) = accumulateState(agg.(s.label), s);
    end

    agg.switching = finalizeState(agg.switching, Parameters.MinStateLengthS);
    agg.stable    = finalizeState(agg.stable, Parameters.MinStateLengthS);
end

function a = emptyStateAgg()
    a.lengthS   = 0;
    a.durationS = 0;
    a.sumGxDur  = 0;
    a.sumStatxDur = 0;
    a.present   = false;
    a.avgG      = NaN;
    a.statVal   = NaN;
    a.maxG      = NaN;
    a.minG      = NaN;
    a.longestSegLengthS = 0;   % internal bookkeeping: lengthS of the longest
                               % segment seen so far, used to pick out that
                               % segment's maxG/minG (not exported directly)
end

function a = accumulateState(a, seg)
    a.lengthS     = a.lengthS + seg.lengthS;
    a.durationS   = a.durationS + seg.durationS;
    a.sumGxDur    = a.sumGxDur + seg.avgG * seg.durationS;
    a.sumStatxDur = a.sumStatxDur + seg.statVal * seg.durationS;
    if seg.lengthS > a.longestSegLengthS
        a.longestSegLengthS = seg.lengthS;
        a.maxG = seg.maxG;
        a.minG = seg.minG;
    end
end

function a = finalizeState(a, minLengthS)
    a.present = a.lengthS >= minLengthS && a.durationS > 0;
    if a.durationS > 0
        a.avgG    = a.sumGxDur / a.durationS;
        a.statVal = a.sumStatxDur / a.durationS;
    end
end

%% ========================================================================
function data = ensureTraceAnalyzed(data, idx)
%ENSURETRACEANALYZED Populate data.cache{idx}, data.segCache{idx} and
% data.aggCache{idx} if they haven't been computed yet.

    if isempty(data.cache{idx})
        [d, t, y, g, ok] = readBreakingTrace(data.fileList{idx}, data.Parameters);
        c.d = d; c.t = t; c.y = y; c.g = g; c.ok = ok;
        data.cache{idx} = c;
    end

    if isempty(data.segCache{idx})
        c = data.cache{idx};
        if ~c.ok || numel(c.d) < 10
            data.segCache{idx} = struct('label', {}, 'startIdx', {}, 'endIdx', {}, ...
                                         'statVal', {}, 'avgG', {}, 'maxG', {}, 'minG', {}, ...
                                         'durationS', {}, 'lengthS', {});
        else
            data.segCache{idx} = segmentTrace(c.d, c.t, c.y, c.g, data.Parameters);
        end
        data.aggCache{idx} = aggregateStates(data.segCache{idx}, data.Parameters);
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
            exportTwoStatePerSegment(fig);
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

    data = ensureTraceAnalyzed(data, idx);
    guidata(fig, data);

    c    = data.cache{idx};
    segs = data.segCache{idx};
    agg  = data.aggCache{idx};
    ax   = data.ax;
    cla(ax);
    hold(ax, 'on');

    if ~c.ok || isempty(c.d)
        text(ax, 0.5, 0.5, 'Could not read this trace', ...
             'HorizontalAlignment', 'center', 'Units', 'normalized');
    else
        d = c.d; y = c.y; g = c.g;
        ylo = min(y) - 0.5;
        yhi = max(y) + 0.5;
        dmax = max(d);

        % Shade trimmed-out start/end regions
        [dStart, dEnd] = conductanceTrimBounds(d, g, P);
        if dStart > min(d)
            fill(ax, [min(d) dStart dStart min(d)], [ylo ylo yhi yhi], ...
                 [0.85 0.85 0.85], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end
        if dEnd < dmax
            fill(ax, [dEnd dmax dmax dEnd], [ylo ylo yhi yhi], ...
                 [0.85 0.85 0.85], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        end

        % Overlay classified segments
        for k = 1:numel(segs)
            s = d(segs(k).startIdx); e = d(segs(k).endIdx);
            if strcmp(segs(k).label, 'switching')
                col = [1 0.3 0.3];   % red-ish
            else
                col = [0.3 0.55 1]; % blue-ish
            end
            fill(ax, [s e e s], [ylo ylo yhi yhi], col, 'FaceAlpha', 0.30, 'EdgeColor', 'none');
        end

        plot(ax, d, y, 'k-', 'LineWidth', 1);
        ylim(ax, [ylo, yhi]);
    end

    set(ax, 'Color', 'w');
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'log_{10}(G/G_0)');
    grid(ax, 'on');
    box(ax, 'on');
    hold(ax, 'off');

    [~, fname] = fileparts(data.fileList{idx});
    exTag = '';
    if data.excluded(idx)
        exTag = '   \color{red}[EXCLUDED FROM EXPORT]';
    end
    bothTag = '';
    if agg.switching.present && agg.stable.present
        bothTag = '   \color[rgb]{0 0.5 0}[BOTH STATES]';
    end
    title(ax, sprintf('Trace %d / %d   -   %s%s%s', ...
          idx, data.Nfiles, fname, bothTag, exTag), 'Interpreter', 'tex');

    statName = 'Std';
    nSwitch = sum(strcmp({segs.label}, 'switching'));
    nStable = sum(strcmp({segs.label}, 'stable'));
    fprintf('Trace %d (%s):\n', idx, fname);
    if agg.switching.present
        fprintf('   switching: %s=%.4g   avgG=%.4e G0   duration=%.4f s   length=%.4f s   (%d segment(s))\n', ...
                statName, agg.switching.statVal, agg.switching.avgG, agg.switching.durationS, agg.switching.lengthS, nSwitch);
    else
        fprintf('   switching: not present\n');
    end
    if agg.stable.present
        fprintf('   stable:    %s=%.4g   avgG=%.4e G0   duration=%.4f s   length=%.4f s   (%d segment(s))\n', ...
                statName, agg.stable.statVal, agg.stable.avgG, agg.stable.durationS, agg.stable.lengthS, nStable);
    else
        fprintf('   stable:    not present\n');
    end

    nExcl = sum(data.excluded);
    set(fig, 'Name', sprintf(['Two-State Noise Analyzer (Per-Segment)  |  %d excluded  |  ' ...
                         'Arrows: browse | Space: exclude | C: clear | ' ...
                         'G: go to # | E: export | Q: quit'], nExcl));
    drawnow;
end

%% ========================================================================
function exportTwoStatePerSegment(fig)
%EXPORTTWOSTATEPERSEGMENT Analyze every non-excluded trace and write three
% files, with ONE ROW PER INDIVIDUAL SEGMENT (not per aggregated state):
%   1. two_state_summary_per_segment.txt  - one row PER SEGMENT, for
%      traces where BOTH switching and stable are present. The trace-level
%      gate is UNCHANGED from TwoStateNoiseAnalyzer.m: it uses the
%      COMBINED, aggregated length of each state against
%      Parameters.MinStateLengthS, exactly as before -
%      Parameters.MinSegmentDurationOutputS below never affects this gate,
%      it only prunes which individual segments get written out once a
%      trace has already qualified.
%   2. segment_scatter_per_segment.txt    - same rows as (1), the raw data
%      behind the 2D histogram.
%   3. segment_histogram2d_per_segment.txt - a 2D histogram (log-spaced
%      mean-G bins x linear std bins) built from (2)'s data, with one
%      count PER SEGMENT (instead of one count per state per trace as in
%      TwoStateNoiseAnalyzer.m). If Parameters.WeightBinsBySegmentDuration
%      is true, each segment contributes its own Duration_s to its bin
%      instead of a flat count of 1, so longer segments carry
%      proportionally more weight in the histogram.

    data = guidata(fig);
    P = data.Parameters;
    statName    = 'Std';
    statUnit    = 'decades';

    outFolder = fullfile(data.folderPath, 'TwoStateAnalysisPerSegment');
    if ~exist(outFolder, 'dir')
        mkdir(outFolder);
    end

    summaryRows = {};   % filename, traceNum, segIdx, state, std, avgG, durationS, maxG, minG
    scatterRows = {};   % filename, traceNum, segIdx, state, std, avgG, durationS

    traceNum = 0;
    nSkippedRead = 0;
    nSegSkippedShort = 0;

    for idx = 1:data.Nfiles
        if data.excluded(idx)
            continue
        end
        data = ensureTraceAnalyzed(data, idx);
        c    = data.cache{idx};
        segs = data.segCache{idx};
        agg  = data.aggCache{idx};
        [~, fname] = fileparts(data.fileList{idx});

        if ~c.ok
            nSkippedRead = nSkippedRead + 1;
            continue
        end

        % ---- Trace-level gate: UNCHANGED from TwoStateNoiseAnalyzer.m ----
        % Uses the combined/aggregated state lengths vs Parameters.MinStateLengthS.
        % Parameters.MinSegmentDurationOutputS plays NO role here.
        bothPresent = agg.switching.present && agg.stable.present;
        if ~bothPresent
            continue
        end

        traceNum = traceNum + 1;

        % ---- Per-segment output: one row per individual segment, in
        % chronological order, filtered ONLY by Parameters.MinSegmentDurationOutputS ----
        segIdx = 0;
        for k = 1:numel(segs)
            s = segs(k);
            if s.durationS < P.MinSegmentDurationOutputS
                nSegSkippedShort = nSegSkippedShort + 1;
                continue
            end
            segIdx = segIdx + 1;
            stdStr = sprintf('%.6g', s.statVal);
            row = {fname, traceNum, segIdx, s.label, stdStr, ...
                   sprintf('%.6g', s.avgG), sprintf('%.6g', s.durationS), ...
                   sprintf('%.6g', s.maxG), sprintf('%.6g', s.minG)};
            summaryRows(end+1, :) = row; %#ok<AGROW>
            scatterRows(end+1, :) = {fname, traceNum, segIdx, s.label, stdStr, ...
                   sprintf('%.6g', s.avgG), sprintf('%.6g', s.durationS)}; %#ok<AGROW>
        end
    end

    % ---- Write file 1: two_state_summary_per_segment.txt ----
    f1 = fullfile(outFolder, 'two_state_summary_per_segment.txt');
    fid = fopen(f1, 'w');
    fprintf(fid, 'Filename\tTraceNumber\tSegmentIndex\tState\tStd_decades\tAvgG_G0\tDuration_s\tMaxG_G0\tMinG_G0\n');
    for k = 1:size(summaryRows, 1)
        fprintf(fid, '%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n', summaryRows{k, :});
    end
    fclose(fid);

    % ---- Write file 2: segment_scatter_per_segment.txt ----
    f2 = fullfile(outFolder, 'segment_scatter_per_segment.txt');
    fid = fopen(f2, 'w');
    fprintf(fid, 'Filename\tTraceNumber\tSegmentIndex\tState\tStd_decades\tAvgG_G0\tDuration_s\n');
    for k = 1:size(scatterRows, 1)
        fprintf(fid, '%s\t%d\t%d\t%s\t%s\t%s\t%s\n', scatterRows{k, :});
    end
    fclose(fid);

    % ---- Build + write file 3: segment_histogram2d_per_segment.txt ----
    f3 = fullfile(outFolder, 'segment_histogram2d_per_segment.txt');
    nRows = size(scatterRows, 1);
    if nRows > 0
        avgG     = cellfun(@(s) str2double(s), scatterRows(:,6));
        statVals = cellfun(@(s) str2double(s), scatterRows(:,5));
        durS     = cellfun(@(s) str2double(s), scatterRows(:,7));
        avgG(avgG <= 0) = P.minlvl;
        logG = log10(avgG);

        gEdges = safeBinEdges(logG, P.HistNBinsG);
        sEdges = safeBinEdges(statVals, P.HistNBinsStat);

        counts = zeros(P.HistNBinsStat, P.HistNBinsG);
        gBinIdx = min(max(discretize(logG, gEdges), 1), P.HistNBinsG);
        sBinIdx = min(max(discretize(statVals, sEdges), 1), P.HistNBinsStat);

        if P.WeightBinsBySegmentDuration
            weights = durS;              % each segment weighted by its own duration [s]
        else
            weights = ones(nRows, 1);    % flat count of 1 per segment
        end
        for k = 1:nRows
            counts(sBinIdx(k), gBinIdx(k)) = counts(sBinIdx(k), gBinIdx(k)) + weights(k);
        end

        fid = fopen(f3, 'w');
        fprintf(fid, '# 2D histogram: rows = %s bins (%s), columns = mean-G bins (log10 G/G0)\n', statName, statUnit);
        fprintf(fid, '# One count PER SEGMENT (not per state per trace).\n');
        if P.WeightBinsBySegmentDuration
            fprintf(fid, '# Bin counts are WEIGHTED by each segment''s Duration_s (longer segments contribute more).\n');
            countFmt     = '%.6g\t';
            countFmtLast = '%.6g\n';
        else
            fprintf(fid, '# Bin counts are UNWEIGHTED (each segment contributes 1, regardless of duration).\n');
            countFmt     = '%d\t';
            countFmtLast = '%d\n';
        end
        fprintf(fid, '# G_bin_edges_log10\t'); fprintf(fid, '%.6g\t', gEdges); fprintf(fid, '\n');
        fprintf(fid, '# Stat_bin_edges\t'); fprintf(fid, '%.6g\t', sEdges); fprintf(fid, '\n');
        fprintf(fid, '# Count matrix follows (one row per Stat bin, low to high; one column per G bin, low to high)\n');
        for r = 1:P.HistNBinsStat
            fprintf(fid, countFmt, counts(r, 1:end-1));
            fprintf(fid, countFmtLast, counts(r, end));
        end
        fclose(fid);
    else
        fid = fopen(f3, 'w');
        fprintf(fid, '# No qualifying segments (from traces with both states present) to histogram.\n');
        fclose(fid);
    end

    guidata(fig, data);

    fprintf('\n=== Export complete ===\n');
    fprintf('Traces with BOTH states present:                     %d\n', traceNum);
    fprintf('Segments exported:                                   %d\n', size(summaryRows, 1));
    fprintf('Segments skipped (< MinSegmentDurationOutputS):      %d\n', nSegSkippedShort);
    fprintf('Traces skipped (unreadable):                         %d\n', nSkippedRead);
    fprintf('Traces excluded manually:                            %d\n', sum(data.excluded));
    fprintf('Wrote:\n  %s\n  %s\n  %s\n', f1, f2, f3);
end
