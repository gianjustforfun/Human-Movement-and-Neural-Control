function [cycles,runs,counts] = readLhoste(filename)
% READLHOSTE Read original Lhoste kinematics and controller gait cycles.
% Keep individual cycles and continuous histories for later calibration.

opts = detectImportOptions(filename);
needed = {'time','condition','gait_state_left','gait_state_right', ...
    'JointPositions_1','JointPositions_2','JointPositions_3','JointPositions_4'};
opts.SelectedVariableNames = needed;
T = readtable(filename,opts);

t = double(T.time(:));
condition = string(T.condition(:));
states = [double(T.gait_state_left(:)),double(T.gait_state_right(:))];
angles = rad2deg([double(T.JointPositions_1(:)), ...
    double(T.JointPositions_2(:)),double(T.JointPositions_3(:)), ...
    double(T.JointPositions_4(:))]);
clear T

% A time reset, gap, or condition change starts a new recording segment.
breaks = [true; condition(2:end)~=condition(1:end-1) | ...
    diff(t)<=0 | diff(t)>0.050 | ~isfinite(diff(t))];
starts = find(breaks);
stops = [starts(2:end)-1;numel(t)];
labels = {'SM_WALKING','transparent_WALKING'};

cycles = struct('condition',{},'leg',{},'segmentIndex',{}, ...
    'segmentFirstRawRow',{},'cycleNumber',{},'HS_s',{},'TO_s',{}, ...
    'nextHS_s',{},'duration_s',{},'stance_s',{},'swing_s',{}, ...
    'toeOffPhase',{}, ...
    'time_s',{},'kinematics',{},'phase_reference_only',{});
segmentIndex = 0;

for seg = 1:numel(starts)
    first = starts(seg);
    last = stops(seg);
    whichCondition = find(strcmp(condition(first),labels),1);
    if isempty(whichCondition)
        continue;
    end
    ts = t(first:last);
    if numel(ts)<2 || any(~isfinite(ts)) || ts(end)-ts(1)<0.50
        continue;
    end
    segmentIndex = segmentIndex+1;

    % Interpolate kinematics only. Event timestamps stay at the original rate.
    t200 = (ts(1):0.005:ts(end))';
    a = interp1(ts,angles(first:last,:),t200,'linear');
    v = nan(size(a));
    v(6:end,:) = (a(6:end,:)-a(1:end-5,:)) ./ ...
        (t200(6:end)-t200(1:end-5));

    for leg = 1:2
        state = states(first:last,leg);
        ix = find(state(1:end-1)~=state(2:end));
        from = state(ix);
        to = state(ix+1);
        eventTime = ts(ix+1);
        HS = eventTime(from==3 & to==0);
        cols = (2*leg-1):(2*leg);
        signal = [a(:,cols),v(:,cols)];

        for j = 1:numel(HS)-1
            hs = HS(j);
            nextHS = HS(j+1);
            inside = eventTime>hs & eventTime<nextHS;
            f = from(inside);
            u = to(inside);
            if ~isequal(f(:)',[0 1 2]) || ~isequal(u(:)',[1 2 3])
                continue;
            end
            changes = eventTime(inside);
            toe = changes(2);       % Controller 1 -> 2, not anatomical toe-off
            duration = nextHS-hs;
            if duration<0.50 || duration>12.0 || ...
                    toe-hs<0.10 || nextHS-toe<0.10
                continue;
            end

            idx = find(t200>=hs & t200<nextHS);
            if numel(idx)<20
                continue;
            end
            tc = t200(idx);
            z = signal(idx,:);
            if tc(1)-hs>0.025 || nextHS-tc(end)>0.025 || ...
                    any(diff(tc)>0.050) || any(~isfinite(z),'all')
                continue;
            end

            one.condition = labels{whichCondition};
            if leg==1
                one.leg = 'L';
            else
                one.leg = 'R';
            end
            one.segmentIndex = segmentIndex;
            one.segmentFirstRawRow = first;
            one.cycleNumber = j;
            one.HS_s = hs;
            one.TO_s = toe;
            one.nextHS_s = nextHS;
            one.duration_s = duration;
            one.stance_s = toe-hs;
            one.swing_s = nextHS-toe;
            one.toeOffPhase = (toe-hs)/duration;
            one.time_s = tc;
            one.kinematics = single(z);
            one.phase_reference_only = single((tc-hs)/duration);
            cycles(end+1,1) = one; %#ok<AGROW>
        end
    end
end

% E06 orders cycles by segment, event time, and leg, not raw time alone.
ordered = cycles([]);
runs = struct('condition',{},'leg',{},'segmentIndex',{}, ...
    'segmentFirstRawRow',{},'time_s',{},'kinematics',{},'phase_reference_only',{}, ...
    'HS_reference_s',{},'TO_reference_s',{},'cycleNumber',{}, ...
    'sourceCycleRows',{});
counts = zeros(2,1);

for c = 1:2
    C = cycles(strcmp({cycles.condition},labels{c}));
    if isempty(C)
        error('No valid Lhoste cycles for %s in %s.',labels{c},filename);
    end
    sideIndex = strcmp({C.leg},'R')';
    [~,order] = sortrows([[C.segmentIndex]',[C.HS_s]',sideIndex]);
    C = C(order);
    ordered = [ordered;C(:)]; %#ok<AGROW>
    counts(c) = numel(C);

    % Join adjacent cycles of the same leg and original segment only.
    for leg = 1:2
        if leg==1
            side = 'L';
        else
            side = 'R';
        end
        rows = find(strcmp({C.leg},side));
        for k = 1:numel(rows)
            u = rows(k);
            one = C(u);
            tc = one.time_s(:);
            z = one.kinematics;
            p = one.phase_reference_only(:);
            exactPhase = (tc-one.HS_s)/one.duration_s;

            % Repair only single-precision rounding at the cycle endpoint.
            roundedEnd = p==1 & exactPhase>=0 & exactPhase<1 & ...
                abs(double(p)-exactPhase)<=2e-6;
            p(roundedEnd) = single(1)-eps('single')/2;
            if any(~isfinite(p)) || any(p<0 | p>=1) || ...
                    any(abs(double(p)-exactPhase)>2e-6) || ...
                    any(diff(tc)<=0 | diff(tc)>0.0125)
                error('Invalid Lhoste cycle history: %s, segment %d.', ...
                    filename,one.segmentIndex);
            end

            append = false;
            if ~isempty(runs) && strcmp(runs(end).condition,labels{c}) && ...
                    strcmp(runs(end).leg,side) && ...
                    runs(end).segmentIndex==one.segmentIndex
                previous = C(runs(end).sourceCycleRows(end));
                append = abs(previous.nextHS_s-one.HS_s)<=0.015 && ...
                    tc(1)>runs(end).time_s(end) && ...
                    tc(1)-runs(end).time_s(end)<=0.0125;
            end
            if append
                runs(end).time_s = [runs(end).time_s;tc];
                runs(end).kinematics = [runs(end).kinematics;z];
                runs(end).phase_reference_only = ...
                    [runs(end).phase_reference_only;p];
                runs(end).HS_reference_s = [runs(end).HS_reference_s; ...
                    repmat(one.HS_s,numel(tc),1)];
                runs(end).TO_reference_s = [runs(end).TO_reference_s; ...
                    repmat(one.TO_s,numel(tc),1)];
                runs(end).cycleNumber = [runs(end).cycleNumber; ...
                    repmat(one.cycleNumber,numel(tc),1)];
                runs(end).sourceCycleRows(end+1,1) = u;
            else
                r.condition = labels{c};
                r.leg = side;
                r.segmentIndex = one.segmentIndex;
                r.segmentFirstRawRow = one.segmentFirstRawRow;
                r.time_s = tc;
                r.kinematics = z;
                r.phase_reference_only = p;
                r.HS_reference_s = repmat(one.HS_s,numel(tc),1);
                r.TO_reference_s = repmat(one.TO_s,numel(tc),1);
                r.cycleNumber = repmat(one.cycleNumber,numel(tc),1);
                r.sourceCycleRows = u;
                runs(end+1,1) = r; %#ok<AGROW>
            end
        end
    end
end
cycles = ordered;
end
