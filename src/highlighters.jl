
## define highlighters used for factor score tables

# hide when 1.0 for diagonals of matrices
hl_hide1 = Highlighter( 
        (data, i, j) -> (data[i,j] ≈ 1.0),
        crayon"black bold"
    );
# positive
hl_top = Highlighter(
        (data, i, j) -> (data[i,j] >= 0.7 && !(data[i,j] ≈ 1.0)),
        crayon"cyan bold"
    );
hl_mid = Highlighter(
        (data, i, j) -> (data[i,j] < 0.7 && data[i,j] >= 0.4),
        crayon"cyan"
    );
hl_low = Highlighter(
        (data, i, j) -> (data[i,j] < 0.4 && data[i,j] >= 0.3),
        crayon"blue"
    );
# hide low values
hl_hide = Highlighter(
        (data, i, j) -> (data[i,j] < 0.3 && data[i,j] >= -0.3),
        crayon"dark_gray"
    );
# negative
hl_topneg = Highlighter(
        (data, i, j) -> (data[i,j] <= -0.7),
        crayon"red bold"
    );
hl_midneg = Highlighter(
        (data, i, j) -> (data[i,j] > -0.7 && data[i,j] <= -0.4),
        crayon"red"
    );
hl_lowneg = Highlighter(
        (data, i, j) -> (data[i,j] > -0.4 && data[i,j] <= -0.3),
        crayon"magenta"
    );

factor_highlighters = (hl_hide1, hl_top, hl_mid, hl_low, hl_hide, hl_topneg, hl_midneg, hl_lowneg)