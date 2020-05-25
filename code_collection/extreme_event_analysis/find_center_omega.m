function [event_lon, event_lat, tag] = find_center_omega(event_lon, event_lat, omega, N_lon, I)

    [Ny, Nx] = size(omega);
    tag = false;
    i = 1;

    pwd_str = pwd;
    if strfind(pwd_str, 'GFDL')
        I = min(3, I);
    else
        I = min(4, I);
    end

    while i <= I % search in concentric boxes of size (2i+1) by (2i+1) around the center of omega field
        temp_omega = omega((Ny + 1) / 2 - i : (Ny + 1) / 2 + i, (Nx + 1) / 2 - i : (Nx + 1) / 2 + i);
        [indy, indx] = find(temp_omega == min(temp_omega(:))); % find the place of minimum omega value
        indy = indy(1); indx = indx(1); % deal with multiple entries

        % if the minimum is not on the edge, it means that it's a local extreme, break
        if indy > 1 && indy < 2 * i + 1 && indx > 1 && indx < 2 * i + 1
            tag = true;
            break;
        end
        i = i + 1;
    end
    if tag
        event_lon = mod(event_lon + (indx - i - 1) - 1, N_lon) + 1;
        event_lat = event_lat + (indy - i - 1);
        if indx - i - 1 ~= 0 || indy - i - 1 ~= 0
            disp(['event is moved, dlon = ', num2str(indx - i - 1), ', dlat = ', num2str(indy - i - 1)]);
        end
    else
        event_lon = NaN;
        event_lat = NaN;
    end
end
