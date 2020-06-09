function [periodtime, periods, max_periods] = period(x,t)
    periods = zeros(1,length(x(1,:)));
    for i = 1:length(x(1,:))
        [pks, locs] = findpeaks(x(:,i));
        if length(locs) >= 2
            periods(i) = length(locs);
        else
            periods(i) = 0;
        end
    end
    
    max_periods = max(periods);
    
    periodtime = zeros(max_periods,length(x(1,:)));
    
    for i = 1:length(x(1,:))
        [pks, locs] = findpeaks(x(:,i));
        for j = 1:max_periods
            if j < length(locs)
                periodtime(j,i) = (t(locs(j+1))-t(locs(j)));
            else
                periodtime(j,i) = 0;
            end
        end
    end
    
end