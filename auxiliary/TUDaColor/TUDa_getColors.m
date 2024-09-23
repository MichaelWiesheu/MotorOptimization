function codes = TUDa_getColors(cols)
    codes = [];
    for i = 1:numel(cols)
        codes = [codes, TUDa_getColor(cols(i))];
    end
end