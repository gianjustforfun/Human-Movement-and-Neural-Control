function [metrics,waves,counts] = readLhoste(filename,phase)
% READLHOSTE Extract walking cycles from the original Lhoste subject CSV.
% Controller 3->0 marks cycle boundaries and 1->2 marks stance-to-swing.
% These are controller-defined events, not measured anatomical toe-off.

% Import only the eight columns needed for this comparison.
opts = detectImportOptions(filename);
needed = {'time','condition','gait_state_left','gait_state_right', ...
    'JointPositions_1','JointPositions_2','JointPositions_3','JointPositions_4'};
opts.SelectedVariableNames = needed;
T = readtable(filename,opts);

t = double(T.time(:));
condition = string(T.condition(:));
states = [double(T.gait_state_left(:)),double(T.gait_state_right(:))];
angles = (180/pi)*[double(T.JointPositions_1(:)), ...
    double(T.JointPositions_2(:)),double(T.JointPositions_3(:)), ...
    double(T.JointPositions_4(:))];
clear T;

% Original time may restart between trials; process continuous sections.
breaks = [true; condition(2:end)~=condition(1:end-1) | ...
    diff(t)<=0 | diff(t)>0.050 | ~isfinite(diff(t))];
starts = find(breaks);
stops = [starts(2:end)-1;length(t)];

labels = ["SM_WALKING","transparent_WALKING"];
values = {zeros(0,6),zeros(0,6)};
waveSums = zeros(length(phase),4,2);
counts = zeros(2,1);

for seg = 1:length(starts)
    first = starts(seg);
    last = stops(seg);
    whichCondition = find(condition(first)==labels);
    if isempty(whichCondition)
        continue;
    end
    originalTime = t(first:last);
    if length(originalTime)<2 || any(~isfinite(originalTime)) || ...
            originalTime(end)-originalTime(1)<0.5
        continue;
    end

    % Resample joint angles to 200 Hz; keep controller event times original.
    t200 = (originalTime(1):0.005:originalTime(end))';
    a = interp1(originalTime,angles(first:last,:),t200,'linear');
    vel = nan(size(a));
    lag = 5;
    vel(lag+1:end,:) = (a(lag+1:end,:)-a(1:end-lag,:)) ./ ...
        (t200(lag+1:end)-t200(1:end-lag));

    for leg = 1:2
        state = states(first:last,leg);
        ix = find(state(1:end-1)~=state(2:end));
        from = state(ix);
        to = state(ix+1);
        changeTime = originalTime(ix+1);
        HS = changeTime(from==3 & to==0);
        cols = (2*leg-1):(2*leg);
        signal = [a(:,cols),vel(:,cols)];

        for k = 1:(length(HS)-1)
            hs = HS(k);
            nextHS = HS(k+1);
            inside = changeTime>hs & changeTime<nextHS;
            f = from(inside);
            u = to(inside);
            if ~isequal(f(:)',[0 1 2]) || ~isequal(u(:)',[1 2 3])
                continue;
            end
            changes = changeTime(inside);
            toe = changes(2);    % Controller state 1->2
            duration = nextHS-hs;
            stance = toe-hs;
            swing = nextHS-toe;
            if duration<0.5 || duration>12 || stance<0.10 || swing<0.10
                continue;
            end

            idx = find(t200>=hs & t200<nextHS);
            if length(idx)<20
                continue;
            end
            tc = t200(idx);
            x = signal(idx,:);
            if tc(1)-hs>0.025 || nextHS-tc(end)>0.025 || ...
                    any(diff(tc)>0.050) || any(~isfinite(x),'all')
                continue;
            end

            normalizedTime = (tc-hs)/duration;
            c = whichCondition;
            waveSums(:,:,c) = waveSums(:,:,c) + ...
                interp1(normalizedTime,x,phase,'linear','extrap');
            values{c}(end+1,:) = [duration,100*stance/duration, ...
                100*swing/duration,stance/swing, ...
                max(x(:,1))-min(x(:,1)),max(x(:,2))-min(x(:,2))]; 
            counts(c) = counts(c)+1;
        end
    end
end

metrics = nan(2,6);
waves = nan(length(phase),4,2);
for c = 1:2
    if counts(c)==0
        error('Missing valid cycles for %s in %s',labels(c),filename);
    end
    metrics(c,:) = median(values{c},1);
    waves(:,:,c) = waveSums(:,:,c)/counts(c);
end
end
