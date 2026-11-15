function eventAlignedAnalysis(lhosteFolder,resultsFolder)
% EVENTALIGNEDANALYSIS Compare SM and Transparent before controller 1->2.
% Uses original Lhoste CSVs only. The event is controller-defined,
% not an independently measured anatomical toe-off.

%% 1. Settings
files = dir(fullfile(lhosteFolder,'*_ExoGait-Benchmark.csv'));
if isempty(files)
    error('No original Lhoste CSV files found.');
end

labels = {'SM_WALKING','transparent_WALKING'};
shortNames = {'SM','Transparent'};
measureNames = {'HipStart_deg','HipEnd_deg','HipChange_deg', ...
    'HipMeanVelocity_deg_s','HipPeakAbsVelocity_deg_s', ...
    'KneeStart_deg','KneeEnd_deg','KneeChange_deg', ...
    'KneeMeanVelocity_deg_s','KneePeakAbsVelocity_deg_s'};
plotNames = {'Hip angle (deg)','Knee angle (deg)', ...
    'Hip angular velocity (deg/s)','Knee angular velocity (deg/s)'};

% Samples are strictly BEFORE the event. t=0 is shown as a marker.
relativeMs = (-300:5:-5)';
relativeTime = relativeMs/1000;
lag = 0.025;               % Backward difference over 25 ms
nSubjects = length(files);
nPoints = length(relativeMs);

% One waveform and ten summary measures per subject and condition.
subjectCurves = nan(nPoints,4,nSubjects,2);
subjectMetrics = nan(nSubjects,2,10);
subjectRows = cell(nSubjects*2,17);

outputFolder = fullfile(resultsFolder,'event_aligned');
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

fprintf('\nEVENT-ALIGNED KINEMATICS - ORIGINAL LHOSTE DATA\n');
fprintf('Time window: -300 to -5 ms; controller transition at 0 ms.\n');

%% 2. Read each subject and collect valid event windows
for s = 1:nSubjects
    filename = fullfile(lhosteFolder,files(s).name);
    id = char(extractBefore(string(files(s).name),'_'));

    opts = detectImportOptions(filename);
    opts.SelectedVariableNames = {'time','condition', ...
        'gait_state_left','gait_state_right', ...
        'JointPositions_1','JointPositions_2', ...
        'JointPositions_3','JointPositions_4'};
    T = readtable(filename,opts);

    t = double(T.time(:));
    condition = string(T.condition(:));
    states = [double(T.gait_state_left(:)),double(T.gait_state_right(:))];
    angles = (180/pi)*[double(T.JointPositions_1(:)), ...
        double(T.JointPositions_2(:)),double(T.JointPositions_3(:)), ...
        double(T.JointPositions_4(:))];
    clear T

    % Separate recording sections at gaps, time resets and mode changes.
    gaps = diff(t);
    breaks = [true; condition(2:end)~=condition(1:end-1) | ...
        gaps<=0 | gaps>0.050 | ~isfinite(gaps)];
    starts = find(breaks);
    stops = [starts(2:end)-1;length(t)];

    sums = zeros(nPoints,4,2,2);  % time, variable, mode, leg
    measures = zeros(2,2,10);     % mode, leg, summary measure
    counts = zeros(2,2);          % valid windows for each mode/leg
    candidates = zeros(2,2);      % eligible complete controller cycles

    for seg = 1:length(starts)
        first = starts(seg);
        last = stops(seg);
        c = find(strcmp(condition(first),labels));
        if isempty(c)
            continue;
        end

        ts = t(first:last);
        if length(ts)<2 || any(~isfinite(ts)) || ts(end)-ts(1)<0.5
            continue;
        end
        a = angles(first:last,:);

        for leg = 1:2
            state = states(first:last,leg);
            ix = find(state(1:end-1)~=state(2:end));
            from = state(ix);
            to = state(ix+1);
            changeTime = ts(ix+1);
            HS = changeTime(from==3 & to==0);
            cols = (2*leg-1):(2*leg);

            for k = 1:length(HS)-1
                hs = HS(k);
                nextHS = HS(k+1);
                inside = changeTime>hs & changeTime<nextHS;
                f = from(inside);
                u = to(inside);
                if ~isequal(f(:)',[0 1 2]) || ~isequal(u(:)',[1 2 3])
                    continue;
                end

                changes = changeTime(inside);
                eventTime = changes(2);  % Controller state 1 -> 2
                duration = nextHS-hs;
                stance = eventTime-hs;
                swing = nextHS-eventTime;
                if duration<0.5 || duration>12 || stance<0.10 || swing<0.10
                    continue;
                end
                candidates(c,leg) = candidates(c,leg)+1;

                query = eventTime+relativeTime;
                past = query-lag;
                % Include the extra 25 ms required for causal velocities.
                % Keep the entire history inside the current stance.
                if past(1)<hs || past(1)<ts(1) || query(end)>ts(end)
                    continue;
                end

                joint = interp1(ts,a(:,cols),query,'linear');
                jointPast = interp1(ts,a(:,cols),past,'linear');
                velocity = (joint-jointPast)/lag;
                x = [joint(:,1),joint(:,2),velocity(:,1),velocity(:,2)];
                if any(~isfinite(x),'all')
                    continue;
                end

                sums(:,:,c,leg) = sums(:,:,c,leg)+x;
                % Start/end angles, angular change, mean velocity and
                % peak absolute velocity. These are per-event measures.
                v = [joint(1,1),joint(end,1), ...
                    joint(end,1)-joint(1,1),mean(velocity(:,1)), ...
                    max(abs(velocity(:,1))),joint(1,2),joint(end,2), ...
                    joint(end,2)-joint(1,2),mean(velocity(:,2)), ...
                    max(abs(velocity(:,2)))];
                measures(c,leg,:) = measures(c,leg,:)+reshape(v,1,1,10);
                counts(c,leg) = counts(c,leg)+1;
            end
        end
    end

    % Each leg contributes equally, then each subject contributes equally.
    % Stop if a subject does not have both legs in both conditions.
    if any(counts(:)==0)
        error('No valid windows for both legs and conditions in %s.',id);
    end

    for c = 1:2
        leftWave = sums(:,:,c,1)/counts(c,1);
        rightWave = sums(:,:,c,2)/counts(c,2);
        subjectCurves(:,:,s,c) = (leftWave+rightWave)/2;

        leftMetrics = reshape(measures(c,1,:),1,10)/counts(c,1);
        rightMetrics = reshape(measures(c,2,:),1,10)/counts(c,2);
        subjectMetrics(s,c,:) = reshape((leftMetrics+rightMetrics)/2,1,1,10);

        row = 2*(s-1)+c;
        v = reshape(subjectMetrics(s,c,:),1,10);
        subjectRows(row,:) = [{id,labels{c},candidates(c,1), ...
            candidates(c,2),counts(c,1),counts(c,2),sum(counts(c,:))}, ...
            num2cell(v)];
    end

    fprintf('%-5s (%2d/%2d): SM %3d/%3d, TR %3d/%3d windows\n', ...
        id,s,nSubjects,sum(counts(1,:)),sum(candidates(1,:)), ...
        sum(counts(2,:)),sum(candidates(2,:)));
end

%% 3. Save subject measurements and paired SM vs Transparent differences
subjectTable = cell2table(subjectRows,'VariableNames', ...
    [{'Subject','Condition','CandidatesLeft','CandidatesRight', ...
      'WindowsLeft','WindowsRight','TotalWindows'},measureNames]);
writetable(subjectTable,fullfile(outputFolder,'subject_metrics.csv'));

pairedRows = cell(nSubjects*10,5);
groupRows = cell(10,8);
for m = 1:10
    sm = reshape(subjectMetrics(:,1,m),nSubjects,1);
    tr = reshape(subjectMetrics(:,2,m),nSubjects,1);
    d = tr-sm;                  % Paired difference within each subject
    for s = 1:nSubjects
        row = (m-1)*nSubjects+s;
        id = char(extractBefore(string(files(s).name),'_'));
        pairedRows(row,:) = {id,measureNames{m},sm(s),tr(s),d(s)};
    end
    groupRows(m,:) = {measureNames{m},mean(sm),mean(tr), ...
        mean(d),std(d),sum(d>0),sum(d<0),nSubjects};
end
pairedTable = cell2table(pairedRows,'VariableNames', ...
    {'Subject','Measure','SM','Transparent','Difference_TR_minus_SM'});
writetable(pairedTable,fullfile(outputFolder,'paired_differences.csv'));

groupTable = cell2table(groupRows,'VariableNames', ...
    {'Measure','MeanSM','MeanTransparent','MeanDifference', ...
     'SDDifference','PositiveSubjects','NegativeSubjects','NSubjects'});
writetable(groupTable,fullfile(outputFolder,'group_metrics.csv'));

%% 4. Figure 1: event-aligned subject mean +/- between-subject SD
colors = lines(2);
figure('Visible','off','Color','w','Position',[100 100 1100 700]);
for j = 1:4
    subplot(2,2,j); hold on;
    curves = zeros(nPoints,nSubjects,2);
    p = gobjects(2,1);
    for c = 1:2
        curves(:,:,c) = reshape(subjectCurves(:,j,:,c),nPoints,nSubjects);
        avg = mean(curves(:,:,c),2);
        sd = std(curves(:,:,c),0,2);
        xx = [relativeMs;flipud(relativeMs)];
        yy = [avg-sd;flipud(avg+sd)];
        fill(xx,yy,colors(c,:),'FaceAlpha',0.13,'EdgeColor','none', ...
            'HandleVisibility','off');
        p(c) = plot(relativeMs,avg,'Color',colors(c,:), ...
            'LineWidth',2);
    end
    xline(0,'--k','HandleVisibility','off');
    xlim([-300 0]);
    xlabel('Time before controller 1->2 (ms)');
    ylabel(plotNames{j}); grid on;
    if j==1
        legend(p,shortNames,'Location','best');
    end
end
sgtitle('Kinematics before stance-to-swing (mean +/- 1 SD across subjects)');
saveas(gcf,fullfile(outputFolder,'01_event_aligned_kinematics.png'));
close(gcf);

%% 5. Figure 2: subject-paired differences for four simple measures
selected = [3 8 5 10];
figure('Visible','off','Color','w','Position',[100 100 1000 650]);
for j = 1:4
    m = selected(j);
    sm = reshape(subjectMetrics(:,1,m),nSubjects,1);
    tr = reshape(subjectMetrics(:,2,m),nSubjects,1);
    d = tr-sm;
    subplot(2,2,j); hold on;
    scatter(1:nSubjects,d,30,'filled');
    yline(0,'--k');
    yline(mean(d),'-r','LineWidth',1.5);
    xlim([0 nSubjects+1]);
    xlabel('Lhoste subject');
    ylabel('Transparent - SM');
    title(strrep(measureNames{m},'_',' '));
    grid on;
end
sgtitle('Paired subject differences before controller 1->2');
saveas(gcf,fullfile(outputFolder,'02_paired_differences.png'));
close(gcf);

fprintf('FINISHED - EVENT-ALIGNED ANALYSIS\n');
fprintf('%d matched Lhoste subjects; figures and CSVs saved in:\n%s\n', ...
    nSubjects,outputFolder);
end
