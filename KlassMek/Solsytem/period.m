function [periodtime, periods] = period(x,t)
    periodtime = length(x(1,:));
    periods = periodtime;
    for i = 1:length(x(1,:))
        [pks, locs] = findpeaks(x(:,i));
        if length(locs) >= 2
            periodtime(i) = (t(locs(2))-t(locs(1)));
            periods(i) = length(locs);
        else
            periodtime(i) = NaN;
            periods(i) = NaN;
        end
    end
end