function [events_sta, num_event] = remove_by_omega(events_sta, num_event, GT_ZERO, plot_level)

    % remove events that satisfy the following criteria:
    % 1) the column-averaged absolute value of omega_QG is larger than omega_threshold;
    % 2) does not have values at plot_level;
    % 3) if GT_ZERO, omega_QG >= 0 at plot_level;
    N = 0;
    omega_threshold = 10.0;
    for j = 1 : length(events_sta.events_sta(:))
        
        ind = true(size(events_sta.events_sta{j}));
        for n = 1 : length(events_sta.events_sta{j})
            level = events_sta.events_sta{j}(n).event_level;
            if abs(nanmean(events_sta.events_sta{j}(n).omega_QG(:))) > omega_threshold || ...
                    isempty(find(level == plot_level)) || ...
                    (exist('GT_ZERO') && GT_ZERO && events_sta.events_sta{j}(n).omega_QG(level == plot_level) >= 0.0)
                ind(n) = false;
                disp(['j = ', num2str(j), ', n = ', num2str(n), ' was removed: either blowed up or omega_500hPa>0']);
                N = N + 1;
            end
        end

        events_sta.events_sta{j} = events_sta.events_sta{j}(ind);
        num_event.num_event(j) = sum(double(ind));
    
    end

    disp(['A total of ', num2str(N), ' events were/was removed.']);

end

