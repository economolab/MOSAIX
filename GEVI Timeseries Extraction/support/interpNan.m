function x = interpNan(x)
    if (size(x,2)==1)
        x = interpNaN1D(x);
    else
        for i = 1:size(x,2)
            x(:,i) = interpNaN1D(x(:,i));
        end
    end
end

function x= interpNaN1D(x)

    nan_idx_x = find(isnan(x));
    non_nan_idx_x = find(~isnan(x));
    x(nan_idx_x) = interp1(non_nan_idx_x, x(non_nan_idx_x), nan_idx_x, 'linear', 'extrap');
end