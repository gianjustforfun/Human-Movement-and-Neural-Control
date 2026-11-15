%% Biomechanical comparison of walking conditions
% Camargo natural walking vs Lhoste SM and Transparent walking.
% Run this file only. The two functions read the ORIGINAL datasets.
% Each subject contributes equally to the group averages.

clear; close all; clc;

%% 1. Folders and settings
folder = fileparts(mfilename('fullpath'));
project = fileparts(folder);
camargoFolder = fullfile(project,'DATASET - CAMARGO (2021)');
lhosteFolder = fullfile(project,'DATASET - LHOSTE (2026)');
resultsFolder = fullfile(folder,'results');
cacheFolder = fullfile(resultsFolder,'cache');

end
if ~isfolder(cacheFolder)
    mkdir(cacheFolder);
end

redo = false;                   % Set true to read the original data again
phase = (0:100)'/100;           % 101 points over one gait cycle
names = {'Camargo','SM','Transparent'};
colors = lines(3);

subjectValues = cell(3,1);     % One row per subject, six scalar measures
subjectCurves = cell(3,1);     % One mean waveform per subject
meanCurves = nan(101,3,4);    % Phase x condition x joint variable
cycleCounts = zeros(3,1);
subjectRows = cell(0,9);

fprintf('\nEXPERIMENT F - BIOMECHANICAL WALKING COMPARISON\n');
fprintf('Inputs: ORIGINAL Camargo and Lhoste datasets only.\n');

%% 2. Camargo: read original marker events and inverse kinematics
subjects = dir(fullfile(camargoFolder,'AB*'));
subjects = subjects([subjects.isdir]);
keep = false(size(subjects));
for s = 1:length(subjects)
    keep(s) = ~isempty(regexp(subjects(s).name,'^AB\d+$','once'));
end
subjects = subjects(keep);
if isempty(subjects)
    error('No Camargo subject folders found.');
end

subjectValues{1} = nan(length(subjects),6);
subjectCurves{1} = nan(101,4,length(subjects));
for s = 1:length(subjects)
    id = subjects(s).name;
    cacheFile = fullfile(cacheFolder,['camargo_' id '.mat']);

    if ~redo && isfile(cacheFile)
        load(cacheFile,'metrics','wave','nCycles');
    else
        [metrics,wave,nCycles] = readCamargo(camargoFolder,id,phase);
        save(cacheFile,'metrics','wave','nCycles');
    end

    subjectValues{1}(s,:) = metrics;
    subjectCurves{1}(:,:,s) = wave;
    subjectRows(end+1,:) = {names{1},id,nCycles,metrics(1), ...
        metrics(2),metrics(3),metrics(4),metrics(5),metrics(6)}; %#ok<SAGROW>
    cycleCounts(1) = cycleCounts(1) + nCycles;
    fprintf('Camargo %-5s (%2d/%2d): %d cycles\n', ...
        id,s,length(subjects),nCycles);
end
fprintf('Camargo total: %d subjects, %d cycles\n\n', ...
    length(subjects),cycleCounts(1));

%% 3. Lhoste: read original CSV and controller transitions
files = dir(fullfile(lhosteFolder,'*_ExoGait-Benchmark.csv'));
if isempty(files)
    error('No original Lhoste CSV files found.');
end

for g = 2:3
    subjectValues{g} = nan(length(files),6);
    subjectCurves{g} = nan(101,4,length(files));
end

for s = 1:length(files)
    id = char(extractBefore(string(files(s).name),'_'));
    cacheFile = fullfile(cacheFolder,['lhoste_' id '.mat']);

    if ~redo && isfile(cacheFile)
        load(cacheFile,'metrics','waves','counts');
    else
        [metrics,waves,counts] = readLhoste( ...
            fullfile(lhosteFolder,files(s).name),phase);
        save(cacheFile,'metrics','waves','counts');
    end

    for c = 1:2
        g = c+1;               % 2 = SM, 3 = Transparent
        subjectValues{g}(s,:) = metrics(c,:);
        subjectCurves{g}(:,:,s) = waves(:,:,c);
        subjectRows(end+1,:) = {names{g},id,counts(c),metrics(c,1), ...
            metrics(c,2),metrics(c,3),metrics(c,4), ...
            metrics(c,5),metrics(c,6)}; %#ok<SAGROW>
        cycleCounts(g) = cycleCounts(g) + counts(c);
    end
    fprintf('Lhoste  %-5s (%2d/%2d): SM %d, Transparent %d cycles\n', ...
        id,s,length(files),counts(1),counts(2));
end

%% 4. Group averages and CSV tables
for g = 1:3
    for j = 1:4
        meanCurves(:,g,j) = mean(subjectCurves{g}(:,j,:),3);
    end
end

subjectSummary = cell2table(subjectRows,'VariableNames', ...
    {'Condition','Subject','NCycles','CycleDuration_s','Stance_pct', ...
     'Swing_pct','StanceSwingRatio','HipROM_deg','KneeROM_deg'});
writetable(subjectSummary,fullfile(resultsFolder,'subject_summary.csv'));

groupRows = cell(3,9);
for g = 1:3
    V = subjectValues{g};
    groupRows(g,:) = {names{g},size(V,1),cycleCounts(g), ...
        mean(V(:,1)),mean(V(:,2)),mean(V(:,3)), ...
        mean(V(:,4)),mean(V(:,5)),mean(V(:,6))};
end
groupSummary = cell2table(groupRows,'VariableNames', ...
    {'Condition','NSubjects','NCycles','CycleDuration_s', ...
     'Stance_pct','Swing_pct','StanceSwingRatio','HipROM_deg','KneeROM_deg'});
writetable(groupSummary,fullfile(resultsFolder,'group_summary.csv'));

%% 5. Figure 1: boxplots of subject-level cycle timing
figure('Visible','off','Color','w','Position',[100 100 1100 370]);
plotNames = {'Cycle duration (s)','Stance (%)','Stance / swing ratio'};
columns = [1 2 4];
for m = 1:3
    subplot(1,3,m); hold on;
    for g = 1:3
        v = subjectValues{g}(:,columns(m));
        boxchart(g*ones(size(v)),v,'BoxFaceColor',colors(g,:), ...
            'BoxWidth',0.5,'CapWidth',0.25);
    end
    xlim([0.5 3.5]);
    xticks(1:3); xticklabels(names);
    ylabel(plotNames{m}); grid on;
end
sgtitle('Cycle timing - subject-level boxplots');
saveas(gcf,fullfile(resultsFolder,'01_cycle_timing.png'));
close(gcf);

%% 6. Figures 2 and 3: joint angles and angular velocities
variables = {'Hip angle (deg)','Knee angle (deg)', ...
             'Hip angular velocity (deg/s)','Knee angular velocity (deg/s)'};
figureNames = {'02_joint_angles.png','03_angular_velocities.png'};
for fig = 1:2
    figure('Visible','off','Color','w','Position',[100 100 1050 400]);
    for j = 1:2
        v = 2*(fig-1)+j;
        subplot(1,2,j); hold on;
        for g = 1:3
            plot(100*phase,meanCurves(:,g,v), ...
                'LineWidth',2,'Color',colors(g,:));
        end
        xlabel('Gait cycle (%)'); ylabel(variables{v});
        xlim([0 100]); grid on;
        legend(names,'Location','best');
    end
    saveas(gcf,fullfile(resultsFolder,figureNames{fig}));
    close(gcf);
end

%% 7. Figure 4: Transparent minus SM (difference of group means)
figure('Visible','off','Color','w','Position',[100 100 1050 700]);
for v = 1:4
    subplot(2,2,v); hold on;
    difference = meanCurves(:,3,v) - meanCurves(:,2,v);
    plot(100*phase,difference,'LineWidth',2);
    yline(0,'--k');
    xlabel('Gait cycle (%)'); ylabel(variables{v});
    title(variables{v});
    xlim([0 100]); grid on;
end
sgtitle('Difference between group means (Transparent - SM)');
saveas(gcf,fullfile(resultsFolder,'04_transparent_minus_SM.png'));
close(gcf);

%% 8. Kinematics immediately before the controller 1->2 transition
% Separate within-subject SM vs Transparent comparison (original CSVs).
eventAlignedAnalysis(lhosteFolder,resultsFolder);

fprintf('\nFINISHED - EXPERIMENT F\n');
fprintf('Camargo: %d cycles | SM: %d | Transparent: %d\n',cycleCounts);
fprintf('Saved 2 CSV files and 4 figures in:\n%s\n',resultsFolder);
