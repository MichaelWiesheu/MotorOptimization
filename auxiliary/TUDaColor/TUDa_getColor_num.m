function code = TUDa_getColor_num(number)
    number =  mod(number-1, 44) + 1;
    digit = mod(number-1, 11) + 1;
    letters = ["a", "b", "c", "d"];
    letter = letters(ceil(number/11));

    code = TUDa_getColor(num2str(digit)+letter);
end
