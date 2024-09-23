% Simple progressbar
% adjusted from https://github.com/ayyu/sholl-analysis-matlab/blob/master/textprogressbar.m

function progressbar(value)

    persistent progressbarStatus;   % 0: not started |  1: ongoing  | 2: finalized
    persistent progressbarStartTime;
    persistent previousValue;

    numberDotsMax = 10;
    numberPercentageLength = 10;

    if nargin==0
        value = '';
    end

    if isempty(progressbarStatus)
        progressbarStatus = 0;
    end
    if isempty(previousValue)
        previousValue = -1;
    end

    % Update printChar according to input
    if ischar(value)
        printChar = value;
        if progressbarStatus == 0
            printChar = [printChar repmat(' ', 1, numberDotsMax+numberPercentageLength) 9];
        elseif progressbarStatus == 1
            progressbarStatus = 2;
            printChar = [9 printChar ' in ' num2str(toc(progressbarStartTime)) ' seconds \n'];
        end
    elseif isnumeric(value)
        percentageValue = floor(value*100);
        if previousValue == percentageValue
            % Only update progressbar if there is a change in progress
            return
        else
            previousValue = percentageValue;
        end
        percentageChar = [num2str(percentageValue) '%%'];
        percentageChar = [percentageChar repmat(' ', 1, numberPercentageLength-length(percentageChar)-1)];
        nDots = floor(value*numberDotsMax);
        dotOut = ['[' repmat('.', 1, nDots) repmat(' ', 1, numberDotsMax-nDots) ']'];
        printChar = [percentageChar dotOut];
        restoreChar = repmat('\b', 1, length(printChar)-1);
    else
        error('Input must be either a char or a numeric value');
    end

    if progressbarStatus == 0
        fprintf(printChar);
        progressbarStartTime = tic();
        progressbarStatus = 1;
    elseif progressbarStatus == 1
        fprintf(restoreChar);
        fprintf(printChar);
    elseif  progressbarStatus == 2
        fprintf(printChar);
        progressbarStatus = 0;
    end
end