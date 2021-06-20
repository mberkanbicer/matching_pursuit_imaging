function [neighbors] = find_neighbor(ind, degree)
    if (degree < 1)
        disp('Degree should be greater than 1')
    else
        [cols, rows] = meshgrid(-degree : degree, -degree : degree);

        [r, c] = size(cols);

        for m = 1:r
            for n=1:c
                neighbors{m, n} = [rows(m, n) cols(m, n)] + ind;
            end
        end
    end