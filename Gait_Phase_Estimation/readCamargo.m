function [trials,cycleDurations,nCycles] = readCamargo(root,id)
% READCAMARGO Build causal phase windows from the original Camargo files.
% Each trial stays separate for subject-based training and validation.

subjectFolder = fullfile(root,id);
sessions = dir(subjectFolder);
sessions = sessions([sessions.isdir]);

trials = struct('session',{},'name',{},'X',{},'phase',{}, ...
    'time',{},'leg',{},'cycleNumber',{});
cycleDurations = zeros(0,1);
nCycles = 0;

% 200-Hz input, 50-Hz anchors, and one second of past kinematics.
history = -200:4:0;

for d = 1:length(sessions)
    base = fullfile(subjectFolder,sessions(d).name,'levelground');
    leftFolder = fullfile(base,'gcLeft');
    if ~isfolder(leftFolder)
        continue;
    end
    files = dir(fullfile(leftFolder,'*.mat'));

    for f = 1:length(files)
        filename = files(f).name;
        if isempty(regexpi(filename,'^levelground_.*\.mat$','once'))
            continue;
        end

        rightFile = fullfile(base,'gcRight',filename);
        ikFile = fullfile(base,'ik',filename);
        if ~isfile(rightFile) || ~isfile(ikFile)
            error('Missing Camargo marker or IK file: %s / %s',id,filename);
        end

        L = load(fullfile(leftFolder,filename),'data');
        R = load(rightFile,'data');
        markers = {L.data,R.data};

        % Exclude trials with duplicated markers or no valid cycles for either leg.
        if isequaln(L.data.ToeOff,R.data.ToeOff)
            continue;
        end
        events = cell(2,1);
        validTrial = true;
        for leg = 1:2
            m = markers{leg};
            mt = double(m.Header(:));
            hsSignal = double(m.HeelStrike(:));
            toSignal = double(m.ToeOff(:));
            if numel(mt)<3 || numel(hsSignal)~=numel(mt) || ...
                    numel(toSignal)~=numel(mt) || any(~isfinite(mt)) || ...
                    any(diff(mt)<=0)
                error('Invalid Camargo marker times: %s / %s',id,filename);
            end
            HS = mt(find(diff(hsSignal)<-50)+1);
            TO = mt(find(diff(toSignal)<-50)+1);
            events{leg} = zeros(0,3);
            for k = 1:length(HS)-1
                toe = TO(TO>HS(k) & TO<HS(k+1));
                if numel(toe)~=1
                    continue;
                end
                stance = toe-HS(k);
                swing = HS(k+1)-toe;
                if stance>=0.15 && stance<=1.80 && ...
                        swing>=0.10 && swing<=1.50
                    events{leg}(end+1,:) = [HS(k),toe,HS(k+1)]; %#ok<AGROW>
                end
            end
            if isempty(events{leg})
                validTrial = false;
            end
        end
        if ~validTrial
            continue;
        end

        A = load(ikFile,'data');
        ik = A.data;
        t = double(ik.Header(:));
        angles = [double(ik.hip_flexion_l(:)), ...
            double(ik.knee_angle_l(:)),double(ik.hip_flexion_r(:)), ...
            double(ik.knee_angle_r(:))];
        if size(angles,1)~=numel(t) || numel(t)<6 || ...
                any(~isfinite(t)) || any(diff(t)<=0)
            error('Invalid Camargo IK data: %s / %s',id,filename);
        end

        % Estimate angular velocities using a five-sample backward difference.
        vel = nan(size(angles));
        vel(6:end,:) = (angles(6:end,:)-angles(1:end-5,:)) ./ ...
            (t(6:end)-t(1:end-5));

        C = struct('leg',{},'number',{},'hs',{},'toe',{}, ...
            'nextHS',{},'duration',{},'time',{},'kin',{},'phase',{});

        for leg = 1:2
            cols = (2*leg-1):(2*leg);
            signal = [angles(:,cols),vel(:,cols)];
            markersForLeg = events{leg};
            for j = 1:size(markersForLeg,1)
                hs = markersForLeg(j,1);
                toe = markersForLeg(j,2);
                nextHS = markersForLeg(j,3);
                idx = find(t>=hs & t<nextHS);
                if numel(idx)<20
                    continue;
                end
                tc = t(idx);
                z = signal(idx,:);
                if tc(1)-hs>0.025 || nextHS-tc(end)>0.025 || ...
                        any(diff(tc)>0.050) || any(~isfinite(z),'all')
                    continue;
                end

                duration = nextHS-hs;
                one.leg = leg;
                one.number = j;
                one.hs = hs;
                one.toe = toe;
                one.nextHS = nextHS;
                one.duration = duration;
                one.time = tc;
                one.kin = single(z);
                one.phase = single((tc-hs)/duration);
                C(end+1,1) = one; %#ok<AGROW>
                cycleDurations(end+1,1) = duration; %#ok<AGROW>
                nCycles = nCycles+1;
            end
        end

        % Join only adjacent accepted cycles from the same trial and leg.
        Xparts = cell(0,1);
        Pparts = cell(0,1);
        Tparts = cell(0,1);
        Lparts = cell(0,1);
        Nparts = cell(0,1);
        for leg = 1:2
            S = C([C.leg]==leg);
            if isempty(S)
                continue;
            end
            [~,order] = sort([S.hs]);
            S = S(order);

            runs = cell(0,1);
            run = struct('time',[],'kin',[],'phase',[],'number',[]);
            previousEnd = NaN;

            for j = 1:numel(S)
                one = S(j);
                tc = one.time(:);
                z = double(one.kin);
                p = double(one.phase(:));
                valid = numel(tc)>=20 && all(isfinite(z),'all') && ...
                    all(isfinite(p)) && all(diff(tc)>0 & diff(tc)<=0.0125) && ...
                    one.hs<one.toe && one.toe<one.nextHS && ...
                    p(1)>=-0.005 && p(end)<=1.005 && ...
                    max(abs(p-(tc-one.hs)/one.duration))<=0.005;
                if ~valid
                    if ~isempty(run.time)
                        runs{end+1,1} = run; %#ok<AGROW>
                    end
                    run = struct('time',[],'kin',[],'phase',[],'number',[]);
                    previousEnd = NaN;
                    continue;
                end

                connected = ~isempty(run.time) && ...
                    abs(one.hs-previousEnd)<=0.015 && ...
                    tc(1)>run.time(end) && tc(1)-run.time(end)<=0.0125;
                if ~connected && ~isempty(run.time)
                    runs{end+1,1} = run; %#ok<AGROW>
                    run = struct('time',[],'kin',[],'phase',[],'number',[]);
                end
                run.time = [run.time;tc];
                run.kin = [run.kin;z];
                run.phase = [run.phase;p];
                run.number = [run.number;repmat(one.number,numel(tc),1)];
                previousEnd = one.nextHS;
            end
            if ~isempty(run.time)
                runs{end+1,1} = run;
            end

            % Each row contains 51 past/current samples per channel.
            for r = 1:numel(runs)
                run = runs{r};
                if numel(run.time)<201
                    continue;
                end
                anchors = (201:4:numel(run.time))';
                goodTime = abs(run.time(anchors)- ...
                    run.time(anchors-200)-1.000)<=0.025;
                anchors = anchors(goodTime);
                if isempty(anchors)
                    continue;
                end
                n = numel(anchors);
                ix = bsxfun(@plus,anchors,history);
                X = single([ ...
                    reshape(run.kin(ix),n,51), ...
                    reshape(run.kin(ix+numel(run.time)),n,51), ...
                    reshape(run.kin(ix+2*numel(run.time)),n,51), ...
                    reshape(run.kin(ix+3*numel(run.time)),n,51)]);
                phi = single(mod(run.phase(anchors),1));
                phi(phi>=1) = single(0);
                good = all(isfinite(X),2) & isfinite(phi) & phi>=0 & phi<1;
                Xparts{end+1,1} = X(good,:); %#ok<AGROW>
                Pparts{end+1,1} = phi(good); %#ok<AGROW>
                Tparts{end+1,1} = run.time(anchors(good)); %#ok<AGROW>
                Lparts{end+1,1} = repmat(uint8(leg),sum(good),1); %#ok<AGROW>
                Nparts{end+1,1} = run.number(anchors(good)); %#ok<AGROW>
            end
        end

        if isempty(Xparts)
            error('No phase windows in Camargo trial: %s / %s',id,filename);
        end
        [~,trialName] = fileparts(filename);
        oneTrial.session = sessions(d).name;
        oneTrial.name = trialName;
        oneTrial.X = vertcat(Xparts{:});
        oneTrial.phase = vertcat(Pparts{:});
        oneTrial.time = vertcat(Tparts{:});
        oneTrial.leg = vertcat(Lparts{:});
        oneTrial.cycleNumber = vertcat(Nparts{:});
        trials(end+1,1) = oneTrial; %#ok<AGROW>
    end
end

if isempty(trials)
    error('No Camargo phase windows found for %s.',id);
end
end
