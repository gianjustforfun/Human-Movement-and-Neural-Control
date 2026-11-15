function [subjectTable,groupTable,ablationTable] = predictTransitions(cacheFolder,resultsFolder,calibrationSummary,model,redoTransitions)
% PREDICTTRANSITIONS Predict Lhoste controller stance-to-swing transitions.
% Fit personalized corrections on the first K cycles per leg only.
% Later controller events are used for offline scoring, never prediction.

%% 1. Settings and saved results
conditions = {'SM_WALKING','transparent_WALKING'};
subjects = sort(unique(string(calibrationSummary.Subject)));
horizons = [0.100 0.200 0.300];
methods = {'Original','TimeRescaled','PhaseRate', ...
    'PhaseTimingRidge','KinematicRidge'};
Kvalues = [5 10 5 10];
sets = {'post5','post10','post10_common','post10_common'};

subjectFile = fullfile(resultsFolder,'transition_subject_metrics.csv');
groupFile = fullfile(resultsFolder,'transition_group_summary.csv');
ablationFile = fullfile(resultsFolder,'transition_ablation_H200.csv');
qcFile = fullfile(resultsFolder,'transition_training_QC.csv');
if ~redoTransitions && all(cellfun(@isfile, ...
        {subjectFile,groupFile,ablationFile,qcFile}))
    subjectTable = readtable(subjectFile);
    groupTable = readtable(groupFile);
    ablationTable = readtable(ablationFile);
    expectedMethods = string(methods);
    currentMethods = unique(string(subjectTable.Method));
    groupMethods = unique(string(groupTable.Method));
    currentColumns = string(ablationTable.Properties.VariableNames);
    validCache = all(ismember(expectedMethods,currentMethods)) && ...
        all(ismember(expectedMethods,groupMethods)) && ...
        all(ismember(["PhaseTimingRidgeRecall_percent", ...
            "KinematicRidgeRecall_percent"],currentColumns));
    if validCache
        fprintf('\nUsing saved transition results. Set redoTransitions=true to rebuild.\n');
        return;
    end
    fprintf('\nSaved results use outdated labels. Recomputing transition results.\n');
end

fprintf('\nPERSONALIZED CONTROLLER TRANSITION PREDICTION\n');
rows = struct([]);
trainingRows = cell(0,6);

%% 2. Personalized models and matched event predictions
for s = 1:numel(subjects)
    id = char(subjects(s));
    D = load(fullfile(cacheFolder,['lhoste_' id '.mat']),'cycles','runs');

    for c = 1:2
        condition = conditions{c};
        C = D.cycles(strcmp({D.cycles.condition},condition));
        runs = D.runs(strcmp({D.runs.condition},condition));
        phaseFile = fullfile(cacheFolder,'phase_predictions', ...
            ['phase_' id '_' condition '.mat']);
        A = load(phaseFile,'protocols');
        protocols = A.protocols;
        if numel(protocols)~=4 || ~isequal( ...
                protocols(3).predictions.key,protocols(4).predictions.key)
            error('Missing or unmatched Lhoste phase protocols for %s.',id);
        end

        selected = string(calibrationSummary.Subject)==string(id) & ...
            string(calibrationSummary.Condition)==string(condition);
        calibrations = sortrows(calibrationSummary(selected,:),'K');
        if height(calibrations)~=2 || ~isequal(calibrations.K(:),[5;10])
            error('Missing K5/K10 calibration for %s %s.',id,condition);
        end

        % Keep the exact transition phase, not the rounded percentage field.
        left = find(strcmp({C.leg},'L'));
        right = find(strcmp({C.leg},'R'));
        calibrations.TransitionPhaseExact = nan(2,1);
        for ki = 1:2
            K = calibrations.K(ki);
            chosen = C([left(1:K);right(1:K)]);
            calibrations.TransitionPhaseExact(ki) = median([chosen.toeOffPhase]);
        end

        % Fit a separate 3-feature and 7-feature Ridge for each K.
        fits = cell(1,2);
        predictions = cell(1,2);
        for ki = 1:2
            cal = calibrations(ki,:);
            fits{ki} = fitCalibration(runs,C,cal,model);
            trainingRows(end+1,:) = {id,condition,cal.K, ...
                fits{ki}.nRows,fits{ki}.nEvents,fits{ki}.nLegs}; %#ok<AGROW>
            predictions{ki} = predictTime(protocols(ki).predictions, ...
                runs,cal,fits{ki});
        end

        for p = 1:4
            P = protocols(p).predictions;
            ki = 1 + (p==2 || p==4);
            cal = calibrations(ki,:);
            R = predictions{ki};
            if p>=3
                % Rates are computed before selecting common post-10 anchors.
                [ok,indices] = ismember(P.key,R.key);
                if ~all(ok) || numel(unique(indices))~=numel(indices)
                    error('Common anchor alignment failed for %s %s.',id,condition);
                end
                fields = fieldnames(R);
                for f = 1:numel(fields)
                    R.(fields{f}) = R.(fields{f})(indices,:);
                end
            end
            if ~isequal(P.key,R.key)
                error('Prediction anchors do not match phase anchors.');
            end

            phase0 = mod(double(P.originalPhase(:)),1);
            phase1 = mod(double(P.timeRescaledPhase(:)),1);
            threshold = double(cal.TransitionPhaseExact);
            duration = double(cal.MedianDuration_s);
            remaining = {mod(threshold-phase0,1)*duration, ...
                mod(threshold-phase1,1)*duration, ...
                R.ttePhaseRate,R.ttePhaseTimingRidge,R.tteKinematicRidge};
            refractory = max(0.300,0.50*duration);
            trueTTE = double(P.TO_reference_s(:))-double(P.time_s(:));
            near = trueTTE>=0.020 & trueTTE<=0.300;

            for h = 1:numel(horizons)
                for m = 1:numel(methods)
                    d = remaining{m};
                    scores = scoreEvents(P,d,horizons(h),refractory);
                    if any(near)
                        nearError = 1000*median(abs(d(near)-trueTTE(near)));
                    else
                        nearError = NaN;
                    end
                    row = struct('Subject',string(id), ...
                        'Condition',string(condition),'CyclesPerLeg',cal.K, ...
                        'EvaluationSet',string(sets{p}), ...
                        'Method',string(methods{m}), ...
                        'Horizon_ms',round(1000*horizons(h)), ...
                        'NAnchors',numel(P.key), ...
                        'NReferenceEvents',scores.NReferenceEvents, ...
                        'NAlarms',scores.NAlarms, ...
                        'TrueAlarms',scores.TrueAlarms, ...
                        'FalseAlarms',scores.FalseAlarms, ...
                        'EventRecall_percent',scores.EventRecall_percent, ...
                        'Precision_percent',scores.Precision_percent, ...
                        'FalseAlarms_perDeviceMin',scores.FalseAlarms_perDeviceMin, ...
                        'MedianAdvance_ms',scores.MedianAdvance_ms, ...
                        'NearEventMedianAbsTTE_ms',nearError, ...
                        'NPreEventAnchors',sum(near), ...
                        'DeviceExposure_min',scores.DeviceExposure_min);
                    if isempty(rows)
                        rows = row;
                    else
                        rows(end+1,1) = row; %#ok<AGROW>
                    end
                end
            end
        end
        fprintf('%-5s %-20s | training rows K5 %4d, K10 %4d\n', ...
            id,condition,fits{1}.nRows,fits{2}.nRows);
    end
end

%% 3. Subject-equal results
subjectTable = struct2table(rows);
writetable(subjectTable,subjectFile);
trainingTable = cell2table(trainingRows,'VariableNames', ...
    {'Subject','Condition','CyclesPerLeg','TrainingRows', ...
    'TrainingEvents','TrainingLegs'});
writetable(trainingTable,qcFile);

groupRows = struct([]);
for c = 1:2
    for p = 1:4
        for h = 1:numel(horizons)
            for m = 1:numel(methods)
                use = subjectTable.Condition==string(conditions{c}) & ...
                    subjectTable.CyclesPerLeg==Kvalues(p) & ...
                    subjectTable.EvaluationSet==string(sets{p}) & ...
                    subjectTable.Horizon_ms==round(1000*horizons(h)) & ...
                    subjectTable.Method==string(methods{m});
                T = subjectTable(use,:);
                if height(T)~=numel(subjects)
                    error('Expected one event result per subject.');
                end
                g = struct('Condition',string(conditions{c}), ...
                    'CyclesPerLeg',Kvalues(p),'EvaluationSet',string(sets{p}), ...
                    'Horizon_ms',round(1000*horizons(h)), ...
                    'Method',string(methods{m}),'NSubjects',height(T), ...
                    'TotalReferenceEvents',sum(T.NReferenceEvents), ...
                    'MeanRecall_percent',mean(T.EventRecall_percent,'omitnan'), ...
                    'MeanPrecision_percent',mean(T.Precision_percent,'omitnan'), ...
                    'MeanFA_perDeviceMin',mean(T.FalseAlarms_perDeviceMin), ...
                    'MedianAdvance_ms',median(T.MedianAdvance_ms,'omitnan'), ...
                    'MedianNearEventAbsTTE_ms', ...
                    median(T.NearEventMedianAbsTTE_ms,'omitnan'));
                if isempty(groupRows)
                    groupRows = g;
                else
                    groupRows(end+1,1) = g; %#ok<AGROW>
                end
            end
        end
    end
end
groupTable = struct2table(groupRows);
writetable(groupTable,groupFile);

% Show the primary matched H200/K5 results for every method.
for c = 1:2
    fprintf('\n%s | K5 common anchors, 200-ms horizon\n',conditions{c});
    for m = 1:numel(methods)
        use = groupTable.Condition==string(conditions{c}) & ...
            groupTable.CyclesPerLeg==5 & ...
            groupTable.EvaluationSet=="post10_common" & ...
            groupTable.Horizon_ms==200 & ...
            groupTable.Method==string(methods{m});
        g = groupTable(use,:);
        fprintf('  %-24s recall %5.2f%% | precision %5.2f%% | FA/min %5.2f\n', ...
            methods{m},g.MeanRecall_percent,g.MeanPrecision_percent, ...
            g.MeanFA_perDeviceMin);
    end
end

% Compare the phase-and-timing Ridge and kinematic Ridge on identical common anchors.
ablationRows = struct([]);
for c = 1:2
    for K = [5 10]
        use = subjectTable.Condition==string(conditions{c}) & ...
            subjectTable.CyclesPerLeg==K & ...
            subjectTable.EvaluationSet=="post10_common" & ...
            subjectTable.Horizon_ms==200;
        A = sortrows(subjectTable(use & ...
            subjectTable.Method=="PhaseTimingRidge",:),'Subject');
        B = sortrows(subjectTable(use & ...
            subjectTable.Method=="KinematicRidge",:),'Subject');
        if height(A)~=numel(subjects) || ~isequal(A.Subject,B.Subject) || ...
                ~isequal(A.NAnchors,B.NAnchors) || ...
                ~isequal(A.NReferenceEvents,B.NReferenceEvents)
            error('Phase-timing and kinematic Ridge evaluations are not matched.');
        end
        delta = B.EventRecall_percent-A.EventRecall_percent;
        row = struct('Condition',string(conditions{c}), ...
            'CyclesPerLeg',K,'Horizon_ms',200, ...
            'Subjects',height(A), ...
            'ReferenceEvents',sum(A.NReferenceEvents), ...
            'PhaseTimingRidgeRecall_percent',mean(A.EventRecall_percent), ...
            'KinematicRidgeRecall_percent',mean(B.EventRecall_percent), ...
            'AddedKinematicsDeltaRecall_pp',mean(delta), ...
            'AddedKinematicsDeltaPrecision_pp', ...
            mean(B.Precision_percent-A.Precision_percent), ...
            'AddedKinematicsDeltaFalseAlarms_perMin', ...
            mean(B.FalseAlarms_perDeviceMin-A.FalseAlarms_perDeviceMin), ...
            'HigherRecallSubjects',sum(delta>1e-9), ...
            'LowerRecallSubjects',sum(delta< -1e-9), ...
            'TiedRecallSubjects',sum(abs(delta)<=1e-9));
        if isempty(ablationRows)
            ablationRows = row;
        else
            ablationRows(end+1,1) = row; %#ok<AGROW>
        end
        fprintf('%-20s K=%2d H200 | phase-timing Ridge %.2f%%, kinematic Ridge %.2f%%, ', ...
            conditions{c},K,row.PhaseTimingRidgeRecall_percent, ...
            row.KinematicRidgeRecall_percent);
        fprintf('delta recall %+.2f pp, delta FA/min %+.2f\n', ...
            row.AddedKinematicsDeltaRecall_pp, ...
            row.AddedKinematicsDeltaFalseAlarms_perMin);
    end
end
ablationTable = struct2table(ablationRows);
writetable(ablationTable,ablationFile);

%% 4. Two figures for the final project
names = {'Original','Rescaled','Phase rate','3-feature Ridge','7-feature Ridge'};
recall = nan(5,2);
for c = 1:2
    for m = 1:5
        use = groupTable.Condition==string(conditions{c}) & ...
            groupTable.CyclesPerLeg==5 & ...
            groupTable.EvaluationSet=="post10_common" & ...
            groupTable.Horizon_ms==200 & ...
            groupTable.Method==string(methods{m});
        recall(m,c) = groupTable.MeanRecall_percent(use);
    end
end
f = figure('Visible','off','Color','w');
bar(recall');
set(gca,'XTick',1:2,'XTickLabel',{'SM','Transparent'});
ylabel('Event recall (%)');
legend(names,'Location','bestoutside');
title('Controller transition prediction, 200 ms, K=5');
grid on;
exportgraphics(f,fullfile(resultsFolder,'transition_recall.png'),'Resolution',180);
close(f);

f = figure('Visible','off','Color','w');
changes = ablationTable(ablationTable.CyclesPerLeg==5,:);
subplot(1,2,1);
bar(changes.AddedKinematicsDeltaRecall_pp);
set(gca,'XTick',1:2,'XTickLabel',{'SM','Transparent'});
ylabel('Recall change (percentage points)');
title('Additional kinematic features');
grid on;
subplot(1,2,2);
bar(changes.AddedKinematicsDeltaFalseAlarms_perMin);
set(gca,'XTick',1:2,'XTickLabel',{'SM','Transparent'});
ylabel('False alarms / device-min (change)');
title('Cost in false alarms');
grid on;
exportgraphics(f,fullfile(resultsFolder,'transition_ablation.png'),'Resolution',180);
close(f);

fprintf('Transition prediction complete. Results saved in results/.\n');
end

function F = fitCalibration(runs,cycles,cal,model)
% Build training examples from the first K completed cycles per leg.
K = cal.K;
left = find(strcmp({cycles.leg},'L'));
right = find(strcmp({cycles.leg},'R'));
selected = [left(1:K);right(1:K)];
selectedCycles = cycles(selected);
offsets = (-50:0)*0.020;
XX = zeros(0,7);
YY = zeros(0,1);
allKeys = zeros(0,1);
legs = zeros(0,1);
nCandidates = 0;

for rr = 1:numel(runs)
    R = runs(rr);
    t = double(R.time_s(:));
    Z = double(R.kinematics);
    if numel(t)<=200 || size(Z,2)~=4 || ...
            any(diff(t)<=0 | diff(t)>0.0125)
        continue;
    end
    ii = (201:4:numel(t))';
    allowed = selectedCycles([selectedCycles.segmentIndex]==R.segmentIndex & ...
        strcmp({selectedCycles.leg},R.leg));
    if isempty(allowed)
        continue;
    end
    hs = [allowed.HS_s]';
    ii = ii(ismember(double(R.HS_reference_s(ii)),hs));
    if isempty(ii)
        continue;
    end
    valid = t(ii)-max(1,double(cal.TemporalScale))>=t(1)-1e-9 & ...
        abs(t(ii)-t(ii-200)-1.0)<=0.015;
    ii = ii(valid);
    if isempty(ii)
        continue;
    end

    % Compute the causal phase rate before thinning calibration examples.
    phi = nan(numel(ii),1);
    for start = 1:256:numel(ii)
        pick = start:min(start+255,numel(ii));
        X = warpedFeatures(t,Z,ii(pick),double(cal.TemporalScale),offsets);
        V = (X(:,model.featureIndices)-double(model.featureMean)) ./ ...
            double(model.featureStd);
        C = V*double(model.coefficients)+double(model.meanCircularTarget);
        phi(pick) = mod(atan2(C(:,2),C(:,1))/(2*pi),1);
    end
    rate = causalRate(t(ii),phi,double(cal.MedianDuration_s));
    keep = 1:3:numel(ii);
    ii = ii(keep);
    phi = phi(keep);
    rate = rate(keep);
    distance = mod(double(cal.TransitionPhaseExact)-phi,1);
    a = distance./rate;
    feat = [a,distance,rate,Z(ii,:)];
    tte = double(R.TO_reference_s(ii))-t(ii);
    eligible = tte>=0.040 & tte<=1.00 & all(isfinite(feat),2);
    if any(eligible)
        XX = [XX;feat(eligible,:)]; %#ok<AGROW>
        YY = [YY;tte(eligible)-a(eligible)]; %#ok<AGROW>
        allKeys = [allKeys;double(R.segmentIndex)*1e9 + ...
            double(strcmp(R.leg,'R'))*1e8 + ...
            double(R.HS_reference_s(ii(eligible)))*1e3]; %#ok<AGROW>
        legs = [legs;repmat(double(strcmp(R.leg,'R'))+1,sum(eligible),1)]; %#ok<AGROW>
    end
    nCandidates = nCandidates+numel(ii);
end

if size(XX,1)<25 || numel(unique(allKeys))<4 || numel(unique(legs))<2
    error('Insufficient calibration data for K=%d.',K);
end

% Standardize using calibration samples only; the Ridge target is residual TTE.
F = struct();
F.nRows = size(XX,1);
F.nEvents = numel(unique(allKeys));
F.nLegs = numel(unique(legs));
F.nCandidateAnchors = nCandidates;
F.maxTTEForTraining_s = 1.0;
F.kinematicRidge = fitRidge(XX,YY,2.0);
F.phaseTimingRidge = fitRidge(XX(:,1:3),YY,2.0);
end

function M = fitRidge(X,y,lambda)
% Both models use exactly the same rows, targets and regularization.
M.xMean = mean(X,1);
M.xStd = max(std(X,0,1),1e-6);
Z = (X-M.xMean)./M.xStd;
M.yMean = mean(y);
M.coeff = (Z'*Z+size(Z,1)*lambda*eye(size(Z,2)))\(Z'*(y-M.yMean));
end

function X = warpedFeatures(t,Z,ii,scale,offsets)
% Kinematics at current and past times only; velocities transform with scale.
n = numel(ii);
query = bsxfun(@plus,t(ii),scale*offsets);
values = interp1(t,Z,query(:),'linear',NaN);
if any(~isfinite(values),'all')
    error('A calibration window crossed a recording boundary.');
end
X = zeros(n,204);
for channel = 1:4
    data = reshape(values(:,channel),n,51);
    if channel>=3
        data = scale*data;
    end
    X(:,(channel-1)*51+(1:51)) = data;
end
end

function rate = causalRate(t,phase,T)
% Use only a 160-ms history and the previously fixed rate limits.
t = double(t(:));
phase = double(phase(:));
prior = 1/T;
rate = prior*ones(size(phase));
j = 0;
for i = 2:numel(phase)
    target = t(i)-0.160+1e-9;
    while j+1<i && t(j+1)<=target
        j = j+1;
    end
    if j==0 || t(i)-t(i-1)>0.045
        rate(i) = prior;
        continue;
    end
    change = mod(phase(i)-phase(j)+0.5,1)-0.5;
    raw = change/(t(i)-t(j));
    rate(i) = min(2.5*prior,max(0.4*prior,0.50*prior+0.50*raw));
end
end

function R = predictTime(P,runs,cal,F)
% Apply both personalized models without using future controller labels.
n = numel(P.key);
R = struct('key',P.key,'ttePhaseRate',zeros(n,1), ...
    'ttePhaseTimingRidge',zeros(n,1),'tteKinematicRidge',zeros(n,1));
runIds = unique(double(P.runIndex(:)),'stable');
for r = 1:numel(runIds)
    rr = runIds(r);
    index = find(double(P.runIndex(:))==rr);
    [t,order] = sort(double(P.time_s(index)));
    index = index(order);
    phi = mod(double(P.timeRescaledPhase(index)),1);
    rate = causalRate(t,phi,double(cal.MedianDuration_s));
    distance = mod(double(cal.TransitionPhaseExact)-phi,1);
    a = min(8.0,distance./rate);
    samples = double(P.anchorSample(index));
    if any(samples>numel(runs(rr).time_s)) || ...
            any(abs(t-double(runs(rr).time_s(samples)))>1e-7)
        error('Kinematics are not aligned with phase anchors.');
    end
    kin = double(runs(rr).kinematics(samples,:));
    features = [a,distance,rate,kin];

    z7 = (features-F.kinematicRidge.xMean)./F.kinematicRidge.xStd;
    correction7 = F.kinematicRidge.yMean+z7*F.kinematicRidge.coeff;
    correction7 = max(-0.35,min(0.35,correction7));
    z3 = (features(:,1:3)-F.phaseTimingRidge.xMean)./F.phaseTimingRidge.xStd;
    correction3 = F.phaseTimingRidge.yMean+z3*F.phaseTimingRidge.coeff;
    correction3 = max(-0.35,min(0.35,correction3));

    b = a;
    p = a;
    eligible = a<=F.maxTTEForTraining_s;
    b(eligible) = max(0,min(8.0,a(eligible)+correction7(eligible)));
    p(eligible) = max(0,min(8.0,a(eligible)+correction3(eligible)));
    R.ttePhaseRate(index) = a;
    R.ttePhaseTimingRidge(index) = p;
    R.tteKinematicRidge(index) = b;
end
end

function M = scoreEvents(P,remaining,H,refractory)
% Controller 1->2 timestamps are used only to score causal alarms.
allRun = double(P.runIndex(:));
tAll = double(P.time_s(:));
legAll = double(P.legIndex(:));
segAll = double(P.segmentIndex(:));
cycAll = double(P.cycleNumber(:));
toAll = double(P.TO_reference_s(:));
nAlarms = 0;
nTrue = 0;
nFalse = 0;
nEvents = 0;
leads = [];

for rr = unique(allRun,'stable')'
    ind = find(allRun==rr);
    [t,order] = sort(tAll(ind));
    ind = ind(order);
    if any(diff(t)<=0) || numel(unique(legAll(ind)))~=1 || ...
            numel(unique(segAll(ind)))~=1
        error('Invalid run identity or time order.');
    end
    positive = double(remaining(ind))<=H+1e-12;
    alarmIdx = zeros(0,1);
    consecutive = 0;
    lastAlarm = -inf;
    for i = 1:numel(t)
        if positive(i)
            if i>1 && t(i)-t(i-1)<=0.045 && positive(i-1)
                consecutive = consecutive+1;
            else
                consecutive = 1;
            end
        else
            consecutive = 0;
        end
        if consecutive>=2 && t(i)-lastAlarm>=refractory-1e-10
            alarmIdx(end+1,1) = i; %#ok<AGROW>
            lastAlarm = t(i);
            consecutive = 0;
        end
    end
    alarms = t(alarmIdx);

    events = unique([cycAll(ind),toAll(ind)],'rows','stable');
    eventTimes = zeros(0,1);
    for e = 1:size(events,1)
        toe = events(e,2);
        belongs = cycAll(ind)==events(e,1);
        eligible = any(belongs & t>=toe-H-1e-9 & t<=toe-0.020+1e-9);
        if eligible
            eventTimes(end+1,1) = toe; %#ok<AGROW>
        end
    end
    eventTimes = sort(eventTimes);
    matched = false(size(eventTimes));
    for i = 1:numel(alarms)
        j = find(~matched & eventTimes>=alarms(i)+0.020-1e-9 & ...
            eventTimes<=alarms(i)+H+1e-9,1,'first');
        if ~isempty(j)
            matched(j) = true;
            nTrue = nTrue+1;
            leads(end+1,1) = eventTimes(j)-alarms(i); %#ok<AGROW>
        else
            nFalse = nFalse+1;
        end
    end
    nAlarms = nAlarms+numel(alarms);
    nEvents = nEvents+numel(eventTimes);
end

% Count overlapping left/right anchor bins once for device exposure.
clockBins = unique([segAll,round(tAll/0.020)],'rows');
deviceMin = size(clockBins,1)*0.020/60;
if nEvents>0
    recall = 100*nTrue/nEvents;
else
    recall = NaN;
end
if nAlarms>0
    precision = 100*nTrue/nAlarms;
else
    precision = 0;
end
if ~isempty(leads)
    advance = 1000*median(leads);
else
    advance = NaN;
end
M = struct('NReferenceEvents',nEvents,'NAlarms',nAlarms, ...
    'TrueAlarms',nTrue,'FalseAlarms',nFalse, ...
    'EventRecall_percent',recall,'Precision_percent',precision, ...
    'FalseAlarms_perDeviceMin',nFalse/deviceMin, ...
    'MedianAdvance_ms',advance,'DeviceExposure_min',deviceMin);
end
