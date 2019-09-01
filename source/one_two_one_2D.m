function sv = one_two_one_2D(v, periodic_d)

    % apply a 1-2-1 filter in two dimensions

    v2 = zeros(size(v) + [2, 2]);
    v2(2:end-1, 2:end-1) = v;
    v2([1,end], 2:end-1) = [v(end, :); v(1, :)];
    v2(2:end-1, [1,end]) = [v(:, end), v(:, 1)];

    sv        = zeros(size(v));
    nv        = min(size(v));

    if nv<3
        error('too few points for filter implementation')
    end

    [w_left, w_right, w_up, w_down] = deal(ones(size(v)));
  
    w_left(:, 2:end) = 1;
    w_left([true(size(w_left, 1), 1), isnan(v(:, 2:end))]) = 0;

    w_right(:, 1:end-1) = 1;
    w_right([isnan(v(:, 1:end-1)), true(size(w_right, 1), 1)]) = 0;

    w_up(2:end, :) = 1;
    w_up([true(1, size(w_up, 2)); isnan(v(2:end, :))]) = 0;

    w_down(1:end-1, :) = 1;
    w_down([isnan(v(1:end-1, :)); true(1, size(w_down, 2))]) = 0;

    w_middle = ones(size(v)) * 4;
    
    if periodic_d == 1
        w_up    (1, :)   = 1;
        w_down  (end, :) = 1;
    elseif periodic_d == 2
        w_left  (:, 1)   = 1;
        w_right (:, end) = 1;
    end

    nan_ind = isnan(v);
    v(isnan(v)) = 0;
    v2(isnan(v2)) = 0;
    sv = (v .* w_middle + ...
          v2(2:end-1, 1:end-2).*w_left + v2(2:end-1, 3:end).*w_right + ...
          v2(1:end-2, 2:end-1).*w_up   + v2(3:end, 2:end-1).*w_down) ./ ...
          (w_left + w_right + w_up + w_down + w_middle);

    sv(nan_ind) = NaN;
    %for j=2:nv-1
    %    sv(j) = 0.5*v(j)+0.25*v(j+1)+0.25*v(j-1);
    %end

    %sv(1) = 0.75*v(1)+0.25*v(2);
    %sv(nv) = 0.75*v(nv)+0.25*v(nv-1);

    %if periodic_d == 1
    %    sv = sv(2:end-1, :);
    %elseif periodic_d == 2
    %    sv = sv(:, 2:end-1);
    %end 





