function TUDa_plotColors()
    cols = ["0a", "1a", "2a", "3a", "4a", "5a", "6a", "7a", "8a", "9a", "10a", "11a";
            "0b", "1b", "2b", "3b", "4b", "5b", "6b", "7b", "8b", "9b", "10b", "11b";
            "0c", "1c", "2c", "3c", "4c", "5c", "6c", "7c", "8c", "9c", "10c", "11c";
            "0d", "1d", "2d", "3d", "4d", "5d", "6d", "7d", "8d", "9d", "10d", "11d"];
    figure()
    hold on

    for iNumber = 1:size(cols,2)
        for iLetter = 1:size(cols,1)
            code = TUDa_getColor(cols(iLetter, iNumber));
            rectangle('Position',[iNumber,5-iLetter,1,1], 'FaceColor',code);
            text(iNumber+0.5, 5-iLetter+0.7, 0, cols(iLetter,iNumber),'HorizontalAlignment','center')
            text(iNumber+0.5, 5-iLetter+0.4, 0, code,'HorizontalAlignment','center','FontSize',7)
        end
    end
end