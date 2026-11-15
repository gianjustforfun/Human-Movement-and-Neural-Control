function [subjectTable,groupTable] = trainLhosteModel(cacheFolder,resultsFolder,camargoModel,calibrationSummary,redoLhoste)
% TRAINLHOSTEMODEL Compare target-domain training with Camargo transfer.
% Train on 14 Lhoste subjects and test on five unseen subjects.
% Controller-defined phase is a training label, never an inference input.

%% 1. Fixed settings and subject-disjoint split
conditions = {'SM_WALKING','transparent_WALKING'};
subjects = sort(unique(string(calibrationSummary.Subject)));
if numel(subjects)~=19 || camargoModel.history_ms~=1000 || ...
        camargoModel.inputColumnCount~=204
    error('Unexpected Lhoste subjects or Camargo feature settings.');
end
rng(574,'twister');
testSubjects = subjects(sort(randperm(numel(subjects),5)));
trainingSubjects = subjects(~ismember(subjects,testSubjects));

modelFile = fullfile(resultsFolder,'lhoste_only_model.mat');
subjectFile = fullfile(resultsFolder,'lhoste_only_subject_metrics.csv');
groupFile = fullfile(resultsFolder,'lhoste_only_group_summary.csv');
figureFile = fullfile(resultsFolder,'lhoste_only_comparison.png');
if ~redoLhoste && isfile(modelFile) && isfile(subjectFile) && ...
        isfile(groupFile) && isfile(figureFile)
    subjectTable = readtable(subjectFile);
    groupTable = readtable(groupFile);
    fprintf('\nUsing saved Lhoste-only comparison. Set redoLhoste=true to rebuild.\n');
    return;
end

fprintf('\nEXPERIMENT E - LHOSTE-ONLY PHASE COMPARISON\n');
fprintf('Training subjects: %s\n',strjoin(cellstr(trainingSubjects),', '));
fprintf('Test subjects: %s\n',strjoin(cellstr(testSubjects),', '));

%% 2. Fit the same 204-feature circular Ridge used for Camargo
% Sample up to 180 equally spaced anchors from each continuous run.
% The two walking conditions are combined for one Lhoste-only model.
nFeatures = 204;
historySamples = 200;               % 1000 ms at 200 Hz
sampleOffsets = -200:4:0;            % 51 samples at 50 Hz
maxPerRun = 180;
lambda = double(camargoModel.lambda); % Keep the selected Camargo lambda fixed

n = 0;
sumX = zeros(1,nFeatures);
sumXX = zeros(nFeatures,nFeatures);
sumC = zeros(1,2);
sumXC = zeros(nFeatures,2);

for s = 1:numel(trainingSubjects)
    id = char(trainingSubjects(s));
    D = load(fullfile(cacheFolder,['lhoste_' id '.mat']),'runs');
    subjectSamples = 0;

    for r = 1:numel(D.runs)
        R = D.runs(r);
        t = double(R.time_s(:));
        if numel(t)<=historySamples
            continue;
        end
        idx = (historySamples+1:4:numel(t))';
        idx = idx(abs(t(idx)-t(idx-historySamples)-1.000)<=0.015);
        if isempty(idx)
            continue;
        end
        idx = idx(unique(round(linspace(1,numel(idx),min(numel(idx),maxPerRun)))));
        count = numel(idx);
        Z = double(R.kinematics);
        samples = bsxfun(@plus,idx,sampleOffsets);
        X = zeros(count,nFeatures);
        for channel = 1:4
            X(:,(channel-1)*51+(1:51)) = reshape(Z(samples,channel),count,51);
        end
        phase = double(R.phase_reference_only(idx)); % Offline supervised labels
        if any(~isfinite(X),'all') || any(~isfinite(phase)) || ...
                any(phase<0 | phase>=1)
            error('Invalid Lhoste training windows for %s.',id);
        end
        C = [cos(2*pi*phase),sin(2*pi*phase)];

        n = n+count;
        subjectSamples = subjectSamples+count;
        sumX = sumX+sum(X,1);
        sumXX = sumXX+X'*X;
        sumC = sumC+sum(C,1);
        sumXC = sumXC+X'*C;
    end
    fprintf('Training %-5s (%2d/%2d): %5d windows\n', ...
        id,s,numel(trainingSubjects),subjectSamples);
end
if n==0
    error('No Lhoste training windows were found.');
end

mu = sumX/n;
sigma = sqrt(max(diag(sumXX)'/n-mu.^2,1e-10));
G = (sumXX-n*(mu'*mu))./(sigma'*sigma);
G = (G+G')/2;
B = (sumXC-mu'*sumC)./sigma';
coefficients = (G+n*lambda*eye(nFeatures))\B;
intercept = sumC/n;

lhosteModel.method = 'circular_phase_ridge';
lhosteModel.history_ms = 1000;
lhosteModel.lambda = lambda;
lhosteModel.featureIndices = 1:nFeatures;
lhosteModel.featureMean = mu;
lhosteModel.featureStd = sigma;
lhosteModel.coefficients = coefficients;
lhosteModel.meanCircularTarget = intercept;
lhosteModel.trainingSubjects = cellstr(trainingSubjects);
lhosteModel.testSubjects = cellstr(testSubjects);
lhosteModel.trainingSamples = n;
lhosteModel.maxSamplesPerRun = maxPerRun;
lhosteModel.conditions = conditions;
save(modelFile,'lhosteModel','-v7');
fprintf('Lhoste-only model trained with %d windows.\n',n);

%% 3. Test all three methods on exactly the same post-K5 anchors
rows = cell(0,7);
chunkSize = 320;
for s = 1:numel(testSubjects)
    id = char(testSubjects(s));
    D = load(fullfile(cacheFolder,['lhoste_' id '.mat']),'runs');

    for c = 1:2
        condition = conditions{c};
        runs = D.runs(strcmp({D.runs.condition},condition));
        predictionFile = fullfile(cacheFolder,'phase_predictions', ...
            ['phase_' id '_' condition '.mat']);
        A = load(predictionFile,'protocols');
        if numel(A.protocols)~=4 || ...
                A.protocols(1).calibrationCyclesPerLeg~=5 || ...
                ~strcmp(A.protocols(1).evaluationSet,'post5')
            error('Unexpected phase evaluation protocol for %s %s.',id,condition);
        end
        P = A.protocols(1).predictions;
        prediction = nan(numel(P.key),1);

        % Recover the same past/current kinematics used by E07 Original.
        for r = 1:numel(runs)
            rowIDs = find(double(P.runIndex)==r);
            if isempty(rowIDs)
                continue;
            end
            t = double(runs(r).time_s(:));
            Z = double(runs(r).kinematics);
            for first = 1:chunkSize:numel(rowIDs)
                chosen = rowIDs(first:min(first+chunkSize-1,numel(rowIDs)));
                anchor = double(P.anchorSample(chosen));
                if any(anchor<=historySamples | anchor>numel(t))
                    error('Invalid phase anchor for %s %s.',id,condition);
                end
                indices = bsxfun(@plus,anchor,sampleOffsets);
                X = zeros(numel(chosen),nFeatures);
                for channel = 1:4
                    X(:,(channel-1)*51+(1:51)) = ...
                        reshape(Z(indices,channel),numel(chosen),51);
                end
                X = (X-lhosteModel.featureMean)./lhosteModel.featureStd;
                C = X*lhosteModel.coefficients+lhosteModel.meanCircularTarget;
                prediction(chosen) = mod(atan2(C(:,2),C(:,1))/(2*pi),1);
            end
        end

        if any(~isfinite(prediction)) || numel(unique(P.key))~=numel(P.key)
            error('Missing or duplicate test predictions for %s %s.',id,condition);
        end
        reference = double(P.referencePhase(:));
        methods = {'Camargo Original','Camargo + K5','Lhoste-only'};
        estimates = {double(P.originalPhase(:)), ...
            double(P.timeRescaledPhase(:)),prediction};

        for m = 1:3
            e = abs(mod(estimates{m}-reference+0.5,1)-0.5);
            rows(end+1,:) = {id,condition,methods{m}, ...
                numel(e),100*mean(e),100*mean(e<=0.10),5}; %#ok<AGROW>
        end
        fprintf('Test %-5s %-20s | Original %.2f, Camargo+K5 %.2f, Lhoste-only %.2f pp\n', ...
            id,condition,rows{end-2,5},rows{end-1,5},rows{end,5});
    end
end

%% 4. Equal-weight test-subject metrics and one comparison figure
subjectTable = cell2table(rows,'VariableNames', ...
    {'Subject','Condition','Method','NAnchors','CircularMAE_pp', ...
    'Within10pp_percent','CalibrationCyclesPerLeg'});
writetable(subjectTable,subjectFile);

groupRows = cell(0,6);
for c = 1:2
    for m = 1:3
        selected = strcmp(subjectTable.Condition,conditions{c}) & ...
            strcmp(subjectTable.Method,methods{m});
        T = subjectTable(selected,:);
        if height(T)~=5 || numel(unique(T.Subject))~=5
            error('Unexpected number of Lhoste test subjects.');
        end
        groupRows(end+1,:) = {conditions{c},methods{m},height(T), ...
            sum(T.NAnchors),mean(T.CircularMAE_pp), ...
            mean(T.Within10pp_percent)}; %#ok<AGROW>
    end
end
groupTable = cell2table(groupRows,'VariableNames', ...
    {'Condition','Method','NSubjects','TotalAnchors', ...
    'CircularMAE_pp','Within10pp_percent'});
writetable(groupTable,groupFile);

values = reshape(groupTable.CircularMAE_pp,3,2)';
fig = figure('Visible','off','Color','w');
bar(values);
set(gca,'XTick',1:2,'XTickLabel',{'SM','Transparent'});
ylabel('Circular phase MAE (percentage points)');
legend(methods,'Location','best');
title('Phase estimation on five unseen Lhoste subjects (post-K5)');
grid on;
exportgraphics(fig,figureFile,'Resolution',180);
close(fig);

fprintf('\nLHOSTE-ONLY COMPARISON - FIVE HELD-OUT SUBJECTS\n');
for c = 1:2
    G = groupTable(strcmp(groupTable.Condition,conditions{c}),:);
    fprintf('%-20s | Original %.2f, Camargo+K5 %.2f, Lhoste-only %.2f pp\n', ...
        conditions{c},G.CircularMAE_pp(1),G.CircularMAE_pp(2), ...
        G.CircularMAE_pp(3));
end
end
