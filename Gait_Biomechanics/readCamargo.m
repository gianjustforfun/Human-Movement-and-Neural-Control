function [metrics,wave,nCycles] = readCamargo(root,id,phase)
% READCAMARGO Extract gait cycles directly from the original Camargo data.
% Marker resets identify heel strike and toe-off; IK gives joint angles.
% Output: six subject medians, four average waveforms, and cycle count.

subjectFolder = fullfile(root,id);
sessions = dir(subjectFolder);
sessions = sessions([sessions.isdir]);

values = zeros(0,6);             % Duration, stance%, swing%, ratio, hip ROM, knee ROM
waveSum = zeros(length(phase),4);
nCycles = 0;

for d = 1:length(sessions)
    base = fullfile(subjectFolder,sessions(d).name,'levelground');
    leftFolder = fullfile(base,'gcLeft');
    if ~isfolder(leftFolder)
        continue;
    end
    trials = dir(fullfile(leftFolder,'*.mat'));

    for f = 1:length(trials)
        filename = trials(f).name;
        if isempty(regexpi(filename,'^levelground_.*\.mat$','once'))
            continue;
        end
        rightFile = fullfile(base,'gcRight',filename);
        ikFile = fullfile(base,'ik',filename);
        if ~isfile(rightFile) || ~isfile(ikFile)
            error('Missing marker or IK file: %s / %s',id,filename);
        end

        L = load(fullfile(leftFolder,filename),'data');
        R = load(rightFile,'data');
        if isequaln(L.data.ToeOff,R.data.ToeOff)
            continue;             % Reject trials with duplicated leg markers
        end
        K = load(ikFile,'data');
        ik = K.data;
        t = double(ik.Header(:));
        if any(~isfinite(t)) || any(diff(t)<=0)
            error('Invalid IK time in %s / %s.',id,filename);
        end
        angles = [double(ik.hip_flexion_l(:)),double(ik.knee_angle_l(:)), ...
                  double(ik.hip_flexion_r(:)),double(ik.knee_angle_r(:))];

        % Compute angular velocities using a five-sample backward difference.
        lag = 5;
        vel = nan(size(angles));
        vel(lag+1:end,:) = (angles(lag+1:end,:) - angles(1:end-lag,:)) ./ ...
            (t(lag+1:end)-t(1:end-lag));

        markers = {L.data,R.data};
        for leg = 1:2
            m = markers{leg};
            mt = double(m.Header(:));
            hsSignal = double(m.HeelStrike(:));
            toSignal = double(m.ToeOff(:));
            if length(mt)~=length(hsSignal) || length(mt)~=length(toSignal) || ...
                    any(diff(mt)<=0)
                error('Invalid marker data in %s / %s.',id,filename);
            end

            HS = mt(find(diff(hsSignal)<-50)+1);
            TO = mt(find(diff(toSignal)<-50)+1);
            % The first two columns are hip/knee of the left leg.
            cols = (2*leg-1):(2*leg);
            signal = [angles(:,cols),vel(:,cols)];

            for k = 1:(length(HS)-1)
                hs = HS(k);
                nextHS = HS(k+1);
                toe = TO(TO>hs & TO<nextHS);
                if numel(toe) ~= 1
                    continue;
                end
                stance = toe-hs;
                swing = nextHS-toe;
                if stance<0.15 || stance>1.80 || swing<0.10 || swing>1.50
                    continue;
                end

                idx = find(t>=hs & t<nextHS);
                if length(idx)<20
                    continue;
                end
                tc = t(idx);
                x = signal(idx,:);
                if tc(1)-hs>0.025 || nextHS-tc(end)>0.025 || ...
                        any(diff(tc)>0.050) || any(~isfinite(x),'all')
                    continue;
                end

                duration = nextHS-hs;
                normalizedTime = (tc-hs)/duration;
                waveSum = waveSum + interp1(normalizedTime,x,phase,'linear','extrap');
                values(end+1,:) = [duration,100*stance/duration, ...
                    100*swing/duration,stance/swing, ...
                    max(x(:,1))-min(x(:,1)),max(x(:,2))-min(x(:,2))]; 
                nCycles = nCycles+1;
            end
        end
    end
end

if nCycles==0
    error('No valid Camargo cycles for %s.',id);
end
metrics = median(values,1);
wave = waveSum/nCycles;
end
