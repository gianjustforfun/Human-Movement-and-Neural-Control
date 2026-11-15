function [subjectTable,groupTable] = evaluatePhase(cacheFolder,resultsFolder,calibrationSummary,model,redoPhase)
% EVALUATEPHASE Test the fixed Camargo phase model on Lhoste walking.
% Compare original and time-rescaled causal histories on matched anchors.

%% 1. Settings
conditions = {'SM_WALKING','transparent_WALKING'};
subjects = sort(unique(string(calibrationSummary.Subject)));
offsets = (-50:0)*0.020;          % 51 samples over the past 1000 ms
history = 200;                   % Original kinematics are at 200 Hz
stride = 4;                      % One prediction every 20 ms
chunkSize = 320;

phaseFolder = fullfile(cacheFolder,'phase_predictions');
if ~isfolder(phaseFolder)
    mkdir(phaseFolder);
end

fprintf('\nEXPERIMENT E - LHOSTE PHASE EVALUATION\n');
rows = cell(0,13);

%% 2. Evaluate each subject and walking condition
for s = 1:numel(subjects)
    id = char(subjects(s));
    file = fullfile(cacheFolder,['lhoste_' id '.mat']);
    D = load(file,'runs');

    for c = 1:2
        condition = conditions{c};
        phaseFile = fullfile(phaseFolder,['phase_' id '_' condition '.mat']);

        if ~redoPhase && isfile(phaseFile)
            A = load(phaseFile,'protocols');
            protocols = A.protocols;
        else
            runs = D.runs(strcmp({D.runs.condition},condition));
            select = string(calibrationSummary.Subject)==string(id) & ...
                string(calibrationSummary.Condition)==string(condition);
            C = sortrows(calibrationSummary(select,:),'K');
            if height(C)~=2 || ~isequal(C.K(:),[5;10])
                error('Missing K5/K10 calibration for %s %s.',id,condition);
            end

            % Separate tests: each K uses only its own calibration cycles.
            P5 = evaluateOne(runs,C(1,:),model,offsets,history,stride,chunkSize);
            P10 = evaluateOne(runs,C(2,:),model,offsets,history,stride,chunkSize);

            % K5 and K10 must use exactly the same post-10 anchors.
            after10 = P5.segmentIndex>C.CutoffSegment(2) | ...
                (P5.segmentIndex==C.CutoffSegment(2) & ...
                P5.time_s>C.CutoffTime_s(2)+1e-9);
            indices5 = find(after10);
            [~,j5,j10] = intersect(P5.key(indices5),P10.key,'stable');
            i5 = indices5(j5);
            if isempty(i5)
                error('No common post-10 anchors for %s %s.',id,condition);
            end
            P5common = pickRows(P5,i5);
            P10common = pickRows(P10,j10);
            if ~isequal(P5common.key,P10common.key) || ...
                    any(P5common.referencePhase~=P10common.referencePhase) || ...
                    any(P5common.originalPhase~=P10common.originalPhase)
                error('K5 and K10 anchor alignment failed for %s %s.',id,condition);
            end

            predictionSets = {P5,P10,P5common,P10common};
            labels = {'post5','post10','post10_common','post10_common'};
            Ks = [5 10 5 10];
            protocols = struct('calibrationCyclesPerLeg',{},'evaluationSet',{}, ...
                'timeScale',{},'predictions',{});
            for p = 1:4
                protocols(p).calibrationCyclesPerLeg = Ks(p);
                protocols(p).evaluationSet = labels{p};
                protocols(p).timeScale = C.TemporalScale(C.K==Ks(p));
                protocols(p).predictions = predictionSets{p};
            end
            save(phaseFile,'protocols','-v7');
        end

        % Calculate summary metrics from the saved paired predictions.
        for p = 1:4
            P = protocols(p).predictions;
            for method = 1:2
                if method==1
                    name = 'Original';
                    predicted = double(P.originalPhase);
                    magnitude = double(P.originalMagnitude);
                else
                    name = 'TimeRescaled';
                    predicted = double(P.timeRescaledPhase);
                    magnitude = double(P.timeRescaledMagnitude);
                end
                reference = double(P.referencePhase);
                signedError = mod(predicted-reference+0.5,1)-0.5;
                err = abs(signedError);
                rows(end+1,:) = {id,condition, ...
                    protocols(p).calibrationCyclesPerLeg, ...
                    protocols(p).evaluationSet,name,numel(reference), ...
                    protocols(p).timeScale,100*mean(err), ...
                    100*mean(err<=0.05),100*mean(err<=0.10), ...
                    100*mean(signedError),100*mean(magnitude<0.20), ...
                    numel(unique(P.key))}; %#ok<AGROW>
            end
        end
        fprintf('%-5s %-20s | K5 %5d, K10 %5d, common %5d anchors\n', ...
            id,condition,numel(protocols(1).predictions.key), ...
            numel(protocols(2).predictions.key), ...
            numel(protocols(3).predictions.key));
    end
end

%% 3. Equal-weight subject averages
subjectTable = cell2table(rows,'VariableNames', ...
    {'Subject','Condition','CyclesPerLeg','EvaluationSet','Method', ...
    'NWindows','TimeScale','CircularMAE_pp','Within5pp_percent', ...
    'Within10pp_percent','SignedBias_pp','LowVector_percent','UniqueAnchors'});
writetable(subjectTable,fullfile(resultsFolder,'phase_subject_metrics.csv'));

groupRows = cell(0,11);
for c = 1:2
    for p = 1:4
        K = [5 10 5 10];
        sets = {'post5','post10','post10_common','post10_common'};
        inSet = strcmp(subjectTable.Condition,conditions{c}) & ...
            subjectTable.CyclesPerLeg==K(p) & ...
            strcmp(subjectTable.EvaluationSet,sets{p});
        original = sortrows(subjectTable(inSet & ...
            strcmp(subjectTable.Method,'Original'),:),'Subject');
        rescaled = sortrows(subjectTable(inSet & ...
            strcmp(subjectTable.Method,'TimeRescaled'),:),'Subject');
        if height(original)~=numel(subjects) || ...
                ~isequal(original.Subject,rescaled.Subject)
            error('Subject-level phase metrics are not paired.');
        end
        groupRows(end+1,:) = {conditions{c},K(p),sets{p}, ...
            numel(subjects),sum(original.NWindows), ...
            mean(original.CircularMAE_pp),mean(rescaled.CircularMAE_pp), ...
            mean(rescaled.CircularMAE_pp-original.CircularMAE_pp), ...
            mean(original.Within10pp_percent), ...
            mean(rescaled.Within10pp_percent), ...
            sum(rescaled.CircularMAE_pp<original.CircularMAE_pp)}; %#ok<AGROW>
    end
end

groupTable = cell2table(groupRows,'VariableNames', ...
    {'Condition','CyclesPerLeg','EvaluationSet','NSubjects', ...
    'TotalMatchedWindows','OriginalMAE_pp','TimeRescaledMAE_pp', ...
    'DeltaMAE_pp','OriginalWithin10_percent', ...
    'TimeRescaledWithin10_percent','TimeRescaledBetterSubjects'});
writetable(groupTable,fullfile(resultsFolder,'phase_group_summary.csv'));

for g = 1:height(groupTable)
    fprintf('%-20s %-13s K=%2d | Original %.2f pp, Rescaled %.2f pp\n', ...
        groupTable.Condition{g},groupTable.EvaluationSet{g}, ...
        groupTable.CyclesPerLeg(g),groupTable.OriginalMAE_pp(g), ...
        groupTable.TimeRescaledMAE_pp(g));
end

%% 4. One figure for the primary K5 comparison
selected = strcmp(groupTable.EvaluationSet,'post5');
G = groupTable(selected,:);
fig = figure('Visible','off','Color','w');
bar([G.OriginalMAE_pp G.TimeRescaledMAE_pp]);
set(gca,'XTick',1:2,'XTickLabel',{'SM','Transparent'});
ylabel('Circular phase MAE (percentage points)');
legend({'Original','Time-rescaled'},'Location','best');
title('Phase estimation after 5 calibration cycles per leg');
grid on;
exportgraphics(fig,fullfile(resultsFolder,'phase_calibration.png'),'Resolution',180);
close(fig);
end

function P = evaluateOne(runs,cal,model,offsets,history,stride,chunkSize)
% Apply both versions of the fixed phase model to the same causal anchors.
scale = double(cal.TemporalScale);
P = struct('key',zeros(0,1,'uint64'), ...
    'referencePhase',zeros(0,1,'single'), ...
    'originalPhase',zeros(0,1,'single'), ...
    'timeRescaledPhase',zeros(0,1,'single'), ...
    'originalMagnitude',zeros(0,1,'single'), ...
    'timeRescaledMagnitude',zeros(0,1,'single'), ...
    'time_s',zeros(0,1),'legIndex',zeros(0,1,'uint8'), ...
    'segmentIndex',zeros(0,1,'uint16'), ...
    'runIndex',zeros(0,1,'uint16'), ...
    'anchorSample',zeros(0,1,'uint32'), ...
    'TO_reference_s',zeros(0,1),'HS_reference_s',zeros(0,1), ...
    'cycleNumber',zeros(0,1));

for rr = 1:numel(runs)
    R = runs(rr);
    t = double(R.time_s(:));
    Z = double(R.kinematics);
    if numel(t)<history+1
        continue;
    end
    if size(Z,1)~=numel(t) || size(Z,2)~=4 || ...
            any(~isfinite(Z(:))) || any(diff(t)<=0 | diff(t)>0.0125)
        error('Invalid Lhoste continuous history.');
    end

    ii = (history+1:stride:numel(t))';
    post = R.segmentIndex>cal.CutoffSegment | ...
        (R.segmentIndex==cal.CutoffSegment & ...
        t(ii)>cal.CutoffTime_s+1e-9);
    ii = ii(post);
    ii = ii(abs(t(ii)-t(ii-history)-1.000)<=0.015);
    ii = ii(t(ii)-max(1,scale)>=t(1)-1e-9);
    if isempty(ii)
        continue;
    end

    for start = 1:chunkSize:numel(ii)
        idx = ii(start:min(start+chunkSize-1,numel(ii)));
        n = numel(idx);
        historyIndex = bsxfun(@plus,idx,-history:4:0);
        X0 = zeros(n,204);
        for channel = 1:4
            X0(:,(channel-1)*51+(1:51)) = ...
                reshape(Z(historyIndex,channel),n,51);
        end

        % Sample older or newer times using only samples up to the anchor.
        queries = bsxfun(@plus,t(idx),scale*offsets);
        values = interp1(t,Z,queries(:),'linear',NaN);
        X1 = zeros(n,204);
        for channel = 1:4
            Q = reshape(values(:,channel),n,51);
            if channel>=3
                Q = scale*Q;     % Convert angular velocity to Camargo time
            end
            X1(:,(channel-1)*51+(1:51)) = Q;
        end
        if any(~isfinite(X0),'all') || any(~isfinite(X1),'all')
            error('A phase history crosses a gap or recording boundary.');
        end

        [phi0,mag0] = phasePrediction(model,X0);
        [phi1,mag1] = phasePrediction(model,X1);
        truth = double(R.phase_reference_only(idx));
        if any(~isfinite(truth) | truth<0 | truth>=1)
            error('Invalid offline phase reference.');
        end

        P.key = [P.key; uint64(rr)*uint64(10000000)+uint64(idx)];
        P.referencePhase = [P.referencePhase; single(truth)];
        P.originalPhase = [P.originalPhase; single(phi0)];
        P.timeRescaledPhase = [P.timeRescaledPhase; single(phi1)];
        P.originalMagnitude = [P.originalMagnitude; single(mag0)];
        P.timeRescaledMagnitude = [P.timeRescaledMagnitude; single(mag1)];
        P.time_s = [P.time_s;t(idx)];
        P.legIndex = [P.legIndex;repmat(uint8(strcmp(R.leg,'R')+1),n,1)];
        P.segmentIndex = [P.segmentIndex;repmat(uint16(R.segmentIndex),n,1)];
        P.runIndex = [P.runIndex;repmat(uint16(rr),n,1)];
        P.anchorSample = [P.anchorSample;uint32(idx)];
        P.TO_reference_s = [P.TO_reference_s;double(R.TO_reference_s(idx))];
        P.HS_reference_s = [P.HS_reference_s;double(R.HS_reference_s(idx))];
        P.cycleNumber = [P.cycleNumber;double(R.cycleNumber(idx))];
    end
end
if isempty(P.key) || numel(unique(P.key))~=numel(P.key)
    error('Empty or duplicate Lhoste phase anchors.');
end
end

function [phase,magnitude] = phasePrediction(model,X)
% The Camargo model is never fitted or modified here.
Z = (double(X(:,model.featureIndices))-double(model.featureMean)) ./ ...
    double(model.featureStd);
C = Z*double(model.coefficients)+double(model.meanCircularTarget);
phase = mod(atan2(C(:,2),C(:,1))/(2*pi),1);
magnitude = hypot(C(:,1),C(:,2));
end

function Q = pickRows(P,idx)
% Keep every prediction and reference field aligned to the same anchors.
Q = P;
fields = fieldnames(P);
for j = 1:numel(fields)
    Q.(fields{j}) = P.(fields{j})(idx,:);
end
end
