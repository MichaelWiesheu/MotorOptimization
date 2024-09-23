function mytimer(name, status)

    persistent timerStruct;
    
    if isempty(timerStruct)
        pos = 0;
    else
        pos = find(vertcat(timerStruct.Name) == name);
        if isempty(pos)
            pos = 0;
        end
    end
    
    if pos == 0
        pos = numel(timerStruct) + 1;
        timerStruct(pos).Name = name;
        timerStruct(pos).StartTime = 0;
        timerStruct(pos).StopTime = 0;
        timerStruct(pos).CalculateTime = 0;
    end
    
    switch status
        case "start"
            timerStruct(pos).StartTime = tic();
        case "stop"
            timerStruct(pos).StopTime = tic();
            timerStruct(pos).CalculateTime = timerStruct(pos).CalculateTime + toc(timerStruct(pos).StartTime);
        case "disp"
            disp("Accumulate time for timer " + name + " : " + num2str(timerStruct(pos).CalculateTime) + " seconds");
        case "dispall"
            for iTimer = 1:numel(timerStruct)
                disp("Accumulate time for timer " + timerStruct(iTimer).Name + " : " + num2str(timerStruct(iTimer).CalculateTime) + " seconds");
            end
        otherwise
            disp("not a valid timer status")
    end
end