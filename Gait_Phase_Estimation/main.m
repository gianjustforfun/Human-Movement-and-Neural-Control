%% Personalized gait-phase calibration
% Run this file only. All input data come from the original datasets.
% Prepares the data, trains the phase model, and predicts controller events.

clear; close all; clc;

%% 1. Folders and settings
folder = fileparts(mfilename('fullpath'));
project = fileparts(folder);
camargoFolder = fullfile(project,'DATASET - CAMARGO (2021)');
lhosteFolder = fullfile(project,'DATASET - LHOSTE (2026)');
resultsFolder = fullfile(folder,'results');
cacheFolder = fullfile(resultsFolder,'cache');

if ~isfolder(cacheFolder)
    mkdir(cacheFolder);
end

redo = false;                   % Set true after changing the reading code
redoModel = false;              % Set true to retrain the Camargo model
redoPhase = false;              % Set true to rebuild Lhoste phase predictions
redoTransitions = false;        % Set true to recompute transition results
runLhosteComparison = true;     % Compare training on Lhoste vs Camargo
redoLhoste = false;             % Set true to repeat the Lhoste-only comparison
Kvalues = [5 10];             % Calibration cycles per leg
conditions = {'SM_WALKING','transparent_WALKING'};

fprintf('\n DATA PREPARATION\n');

%% 2. Camargo: individual cycles and causal phase windows
subjects = dir(fullfile(camargoFolder,'AB*'));
subjects = subjects([subjects.isdir]);
keep = false(size(subjects));
for s = 1:numel(subjects)
    keep(s) = ~isempty(regexp(subjects(s).name,'^AB\d+$','once'));
end
subjects = subjects(keep);

camargoDurations = zeros(0,1);
camargoCycles = 0;
camargoWindows = 0;
camargoTrials = 0;

for s = 1:numel(subjects)
    id = subjects(s).name;
    cacheFile = fullfile(cacheFolder,['camargo_' id '.mat']);

    if ~redo && isfile(cacheFile)
        load(cacheFile,'trials','cycleDurations','nCycles');
    else
        [trials,cycleDurations,nCycles] = readCamargo(camargoFolder,id);
        save(cacheFile,'trials','cycleDurations','nCycles','-v7');
    end

    nWindows = 0;
    for j = 1:numel(trials)
        nWindows = nWindows+size(trials(j).X,1);
    end
    camargoDurations = [camargoDurations;cycleDurations]; %#ok<AGROW>
    camargoCycles = camargoCycles+nCycles;
    camargoWindows = camargoWindows+nWindows;
    camargoTrials = camargoTrials+numel(trials);
    fprintf('Camargo %-5s (%2d/%2d): %4d cycles, %6d windows\n', ...
        id,s,numel(subjects),nCycles,nWindows);
    clear trials cycleDurations
end

% E06 uses the median of all Camargo cycles, not subject medians.
camargoPeriod = median(camargoDurations);
fprintf('Camargo: %d subjects, %d trials, %d cycles, %d windows\n', ...
    numel(subjects),camargoTrials,camargoCycles,camargoWindows);
fprintf('Reference cycle duration: %.4f s\n\n',camargoPeriod);

% These counts were obtained with the original E01-E04 processing.
if numel(subjects)~=22 || camargoTrials~=680 || ...
        camargoCycles~=14942 || camargoWindows~=835934
    error('Camargo counts differ from the verified E01-E04 results.');
end

%% 3. Train and validate the Camargo phase model
modelFile = fullfile(resultsFolder,'phase_model.mat');
modelRebuilt = redo || redoModel || ~isfile(modelFile);
if modelRebuilt
    [model,modelValidation,subjectValidation,trainingInfo] = ...
        trainPhaseModel(cacheFolder,{subjects.name});
    save(modelFile,'model','trainingInfo','modelValidation', ...
        'subjectValidation','-v7');
    writetable(modelValidation,fullfile(resultsFolder,'phase_validation_models.csv'));
    writetable(subjectValidation,fullfile(resultsFolder,'phase_validation_subjects.csv'));
else
    load(modelFile,'model','trainingInfo');
    fprintf('Using saved Camargo phase model. Set redoModel=true to retrain.\n');
end
fprintf('Camargo validation: %.2f pp MAE, %.1f%% within 10 pp\n\n', ...
    trainingInfo.validationMAE_pp,trainingInfo.validationWithin10_percent);

%% 4. Lhoste: individual controller cycles and continuous histories
files = dir(fullfile(lhosteFolder,'*_ExoGait-Benchmark.csv'));
if isempty(files)
    error('No original Lhoste subject CSV files found.');
end

lhosteCounts = zeros(2,1);
calibrationRows = cell(0,9);
for s = 1:numel(files)
    id = char(extractBefore(string(files(s).name),'_'));
    cacheFile = fullfile(cacheFolder,['lhoste_' id '.mat']);

    if ~redo && isfile(cacheFile)
        load(cacheFile,'cycles','runs','counts');
    else
        [cycles,runs,counts] = readLhoste( ...
            fullfile(lhosteFolder,files(s).name));
        save(cacheFile,'cycles','runs','counts','-v7');
    end
    lhosteCounts = lhosteCounts+counts;

    % Initial temporal calibration for each condition
    for c = 1:2
        C = cycles(strcmp({cycles.condition},conditions{c}));
        left = find(strcmp({C.leg},'L'));
        right = find(strcmp({C.leg},'R'));

        for K = Kvalues
            if numel(left)<K || numel(right)<K
                error('Not enough %s calibration cycles for %s.', ...
                    conditions{c},id);
            end
            selected = [left(1:K);right(1:K)];
            segment = [C(selected).segmentIndex]';
            cycleEnd = [C(selected).nextHS_s]';
            lastSegment = max(segment);
            cutoffTime = max(cycleEnd(segment==lastSegment));

            duration = [C(selected).duration_s]';
            scale = median(duration)/camargoPeriod;
            transitionPhase = median([C(selected).toeOffPhase]);

            % Count anchors after BOTH legs finish their initial K cycles.
            nPost = 0;
            for r = 1:numel(runs)
                if ~strcmp(runs(r).condition,conditions{c})
                    continue;
                end
                anchors = (201:4:numel(runs(r).time_s))';
                isPost = runs(r).segmentIndex>lastSegment | ...
                    (runs(r).segmentIndex==lastSegment & ...
                    runs(r).time_s(anchors)>cutoffTime+1e-9);
                nPost = nPost+sum(isPost);
            end
            if nPost==0
                error('No post-calibration anchors for %s, %s, K=%d.', ...
                    id,conditions{c},K);
            end

            calibrationRows(end+1,:) = {id,conditions{c},K, ...
                lastSegment,cutoffTime,median(duration),scale, ...
                100*transitionPhase,nPost}; %#ok<AGROW>
        end
    end

    fprintf('Lhoste  %-5s (%2d/%2d): SM %3d, Transparent %3d cycles\n', ...
        id,s,numel(files),counts(1),counts(2));
    clear cycles runs
end

if numel(files)~=19 || ~isequal(lhosteCounts,[1467;1497])
    error('Lhoste counts differ from the verified E03 results.');
end

calibrationSummary = cell2table(calibrationRows,'VariableNames', ...
    {'Subject','Condition','K','CutoffSegment','CutoffTime_s', ...
    'MedianDuration_s','TemporalScale','TransitionPhase_pct', ...
    'PostCalibrationAnchors'});
writetable(calibrationSummary,fullfile(resultsFolder,'calibration_summary.csv'));

fprintf('\nLhoste: %d subjects, %d SM cycles, %d Transparent cycles\n', ...
    numel(files),lhosteCounts(1),lhosteCounts(2));
fprintf('Prepared %d subject-condition-K calibrations.\n', ...
    height(calibrationSummary));
%% 5. Lhoste: evaluate original and time-rescaled phase
[phaseSubjects,phaseGroups] = evaluatePhase(cacheFolder,resultsFolder, ...
    calibrationSummary,model,redo || redoPhase || modelRebuilt); %#ok<ASGLU>

fprintf('Phase evaluation complete. Results saved in results/.\n');

%% 6. Lhoste: personalized time-to-transition and feature ablation
[transitionSubjects,transitionGroups,ablation] = predictTransitions( ...
    cacheFolder,resultsFolder,calibrationSummary,model, ...
    redo || redoModel || redoPhase || redoTransitions || modelRebuilt); %#ok<ASGLU>

%% 7. Optional Lhoste-only training and held-out comparison
if runLhosteComparison
    [lhosteSubjects,lhosteGroups] = trainLhosteModel( ...
        cacheFolder,resultsFolder,model,calibrationSummary, ...
        redo || redoModel || redoPhase || redoLhoste || modelRebuilt); %#ok<ASGLU>
end

fprintf('Complete.\n');

