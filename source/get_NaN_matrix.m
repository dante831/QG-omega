function NaN_matrix = get_NaN_matrix(data, threshold)

    if threshold <= 0
        error('Need threshold larger than 0!')
    end
    % if any entries in data has value larger than threshold, remove this grid point
    temp_ind = sum(abs(data) > threshold, 3) >= 1; 
    NaN_matrix = ones([size(data, 1), size(data, 2)]);
    NaN_matrix(temp_ind) = NaN;

end
