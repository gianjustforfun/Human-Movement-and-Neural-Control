function [model,modelValidation,subjectValidation,trainingInfo] = trainPhaseModel(cacheFolder,subjectIDs)
% TRAINPHASEMODEL Train and validate the Camargo circular phase Ridge model.
% All features come from the original Camargo data cached by readCamargo.

%% 1. Fixed training settings
subjects = sort(string(subjectIDs(:)));
rng(574,'twister');
validationSubjects = subjects(sort(randperm(numel(subjects),5)));
trainingSubjects = subjects(~ismember(subjects,validationSubjects));

histories = [500 1000];
lambdas = [0.0001 0.001 0.01 0.1];
last500 = 26:51;
featureSets = {[last500 last500+51 last500+102 last500+153],1:204};
maxPerTrial = 180;

fprintf('\nEXPERIMENT E - CAMARGO PHASE MODEL\n');
fprintf('Training subjects: %s\n',strjoin(cellstr(trainingSubjects),', '));
fprintf('Validation subjects: %s\n',strjoin(cellstr(validationSubjects),', '));

%% 2. Accumulate training statistics, one subject at a time
% The second group includes all subjects for the final model refit.
stats = cell(2,2);
for h = 1:2
    p = numel(featureSets{h});
    for g = 1:2
        stats{h,g} = struct('n',0,'sumX',zeros(1,p), ...
            'sumXX',zeros(p,p),'sumC',zeros(1,2),'sumXC',zeros(p,2));
    end
end

nTrials = 0;
nWindows = 0;
for s = 1:numel(subjects)
    file = fullfile(cacheFolder,['camargo_' char(subjects(s)) '.mat']);
    S = load(file,'trials');
    trials = S.trials;
    for f = 1:numel(trials)
        n = size(trials(f).X,1);
        if size(trials(f).X,2)~=204 || numel(trials(f).phase)~=n || n==0
            error('Unexpected Camargo feature dimensions: %s.',file);
        end
        nTrials = nTrials+1;
        nWindows = nWindows+n;

        % E05 uses up to 180 equally spaced windows from each trial.
        idx = unique(round(linspace(1,n,min(n,maxPerTrial))));
        X = double(trials(f).X(idx,:));
        phi = mod(double(trials(f).phase(idx)),1);
        phi = phi(:);
        if any(~isfinite(X),'all') || any(~isfinite(phi))
            error('Invalid Camargo features in %s.',file);
        end
        C = [cos(2*pi*phi),sin(2*pi*phi)];

        groups = 2;
        if ~ismember(subjects(s),validationSubjects)
            groups = [1 2];
        end
        for h = 1:2
            F = X(:,featureSets{h});
            sx = sum(F,1);
            sxx = F'*F;
            sc = sum(C,1);
            sxc = F'*C;
            for g = groups
                A = stats{h,g};
                A.n = A.n+size(F,1);
                A.sumX = A.sumX+sx;
                A.sumXX = A.sumXX+sxx;
                A.sumC = A.sumC+sc;
                A.sumXC = A.sumXC+sxc;
                stats{h,g} = A;
            end
        end
    end
    fprintf('Training data %-5s (%2d/%2d): %3d trials\n', ...
        char(subjects(s)),s,numel(subjects),numel(trials));
end

if numel(subjects)~=22 || nTrials~=680 || nWindows~=835934
    error('Camargo training data differ from the verified E04 counts.');
end

%% 3. Fit the eight models using training subjects only
models = struct('history_ms',{},'lambda',{},'columns',{}, ...
    'mu',{},'sigma',{},'intercept',{},'coefficients',{});
for h = 1:2
    A = stats{h,1};
    n = A.n;
    mu = A.sumX/n;
    sigma = sqrt(max(diag(A.sumXX)'/n-mu.^2,1e-10));
    G = (A.sumXX-n*(mu'*mu))./(sigma'*sigma);
    G = (G+G')/2;
    C = (A.sumXC-mu'*A.sumC)./sigma';
    for k = 1:numel(lambdas)
        j = numel(models)+1;
        models(j).history_ms = histories(h);
        models(j).lambda = lambdas(k);
        models(j).columns = featureSets{h};
        models(j).mu = mu;
        models(j).sigma = sigma;
        models(j).intercept = A.sumC/n;
        models(j).coefficients = (G+n*lambdas(k)*eye(size(G)))\C;
    end
end

%% 4. Validate on every window of the five held-out subjects
nVal = numel(validationSubjects);
nModels = numel(models);
subjectCount = zeros(nVal,1);
errorSum = zeros(nVal,nModels);
within5 = zeros(nVal,nModels);
within10 = zeros(nVal,nModels);
uncertain = zeros(nVal,nModels);

for v = 1:nVal
    file = fullfile(cacheFolder,['camargo_' char(validationSubjects(v)) '.mat']);
    S = load(file,'trials');
    for f = 1:numel(S.trials)
        X = double(S.trials(f).X);
        phi = double(S.trials(f).phase(:));
        n = size(X,1);
        subjectCount(v) = subjectCount(v)+n;
        for j = 1:nModels
            M = models(j);
            Z = (X(:,M.columns)-M.mu)./M.sigma;
            C = Z*M.coefficients+M.intercept;
            predicted = mod(atan2(C(:,2),C(:,1))/(2*pi),1);
            signedError = mod(predicted-phi+0.5,1)-0.5;
            e = abs(signedError);
            errorSum(v,j) = errorSum(v,j)+sum(e);
            within5(v,j) = within5(v,j)+sum(e<=0.05);
            within10(v,j) = within10(v,j)+sum(e<=0.10);
            uncertain(v,j) = uncertain(v,j)+sum(hypot(C(:,1),C(:,2))<0.20);
        end
    end
end

if any(subjectCount==0)
    error('A Camargo validation subject has no windows.');
end
subjectMAE = errorSum./subjectCount;
modelMAE = mean(subjectMAE,1);
modelWithin5 = mean(within5./subjectCount,1);
modelWithin10 = mean(within10./subjectCount,1);
modelUncertain = mean(uncertain./subjectCount,1);
[bestMAE,best] = min(modelMAE);
winner = models(best);

fprintf('\nSUBJECT-DISJOINT VALIDATION\n');
for j = 1:nModels
    fprintf('History %4d ms, lambda %.4f: MAE %.2f pp, within 10 pp %.1f%%\n', ...
        models(j).history_ms,models(j).lambda, ...
        100*modelMAE(j),100*modelWithin10(j));
end
fprintf('Selected: %d ms, lambda %.4f\n',winner.history_ms,winner.lambda);
for v = 1:nVal
    fprintf('%s: %d windows, MAE %.2f pp, within 10 pp %.1f%%\n', ...
        char(validationSubjects(v)),subjectCount(v), ...
        100*subjectMAE(v,best),100*within10(v,best)/subjectCount(v));
end

%% 5. Refit the selected model on all 22 Camargo subjects
h = find(histories==winner.history_ms,1);
A = stats{h,2};
n = A.n;
mu = A.sumX/n;
sigma = sqrt(max(diag(A.sumXX)'/n-mu.^2,1e-10));
G = (A.sumXX-n*(mu'*mu))./(sigma'*sigma);
G = (G+G')/2;
C = (A.sumXC-mu'*A.sumC)./sigma';
coefficients = (G+n*winner.lambda*eye(size(G)))\C;

model.method = 'circular_phase_ridge';
model.history_ms = winner.history_ms;
model.lambda = winner.lambda;
model.featureIndices = featureSets{h};
model.featureMean = mu;
model.featureStd = sigma;
model.coefficients = coefficients;
model.meanCircularTarget = A.sumC/n;
model.inputColumnCount = 204;
model.inputSamplingHz = 50;
model.featureChannelOrder = {'hip_angle_deg','knee_angle_deg', ...
    'hip_velocity_deg_s','knee_velocity_deg_s'};
model.trainingSubjects = cellstr(subjects);
model.validationSubjectsForSelection = cellstr(validationSubjects);
model.trainingSamplesFinal = n;
model.usesLhosteForTraining = false;

modelValidation = table([models.history_ms]',[models.lambda]', ...
    100*modelMAE',100*modelWithin5',100*modelWithin10', ...
    100*modelUncertain', ...
    'VariableNames',{'History_ms','Lambda','CircularMAE_pp', ...
    'Within5pp_percent','Within10pp_percent','LowVectorMagnitude_percent'});
subjectValidation = table(validationSubjects,subjectCount, ...
    100*subjectMAE(:,best),100*within5(:,best)./subjectCount, ...
    100*within10(:,best)./subjectCount, ...
    'VariableNames',{'Subject','Windows','CircularMAE_pp', ...
    'Within5pp_percent','Within10pp_percent'});

trainingInfo.trainingSubjects = cellstr(trainingSubjects);
trainingInfo.validationSubjects = cellstr(validationSubjects);
trainingInfo.randomSeed = 574;
trainingInfo.maxSampledWindowsPerTrial = maxPerTrial;
trainingInfo.histories_ms = histories;
trainingInfo.lambdas = lambdas;
trainingInfo.trials = nTrials;
trainingInfo.totalWindows = nWindows;
trainingInfo.validationWindows = sum(subjectCount);
trainingInfo.trainingSamples = stats{2,1}.n;
trainingInfo.finalSamples = n;
trainingInfo.validationMAE_pp = 100*bestMAE;
trainingInfo.validationWithin10_percent = 100*modelWithin10(best);

fprintf('\nCamargo validation: MAE %.2f pp, within 10 pp %.1f%%\n', ...
    trainingInfo.validationMAE_pp,trainingInfo.validationWithin10_percent);
fprintf('Final model refitted with %d Camargo samples.\n',n);
end
