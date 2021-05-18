function [komsular] = komsubul(ind, derece)
    if (derece < 1)
        disp('derece 1''den buyuk olmali')
    else
        [cols, rows] = meshgrid(-derece : derece, -derece : derece);

        [r, c] = size(cols);

        for m = 1:r
            for n=1:c
                komsular{m, n} = [rows(m, n) cols(m, n)] + ind;
            end
        end
    end