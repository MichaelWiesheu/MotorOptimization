% https://www.mathcha.io/editor/mXDrwcoBi8nsGwhkN6zJwuoPlMW0cMzEE4nc16rnon

function [srf, patches, parameters] = GEO_DISstator_2(parameters)
    draw_geometry = true;
    % GAP = 1e-3;
    % GapType = 'free'  TBD usage in JMAG?
    % drive standard parameters
    HEIGHT = 30e-3;
    POLES = 4;
    SLOTS = 24;
    % stator standard parameters
    RD1 = 100e-3;
    SD1 = 204e-3;
    SD2 = 102e-3;
    ST = 1.6e-3;
    SW1 = 6.0e-3;
    SW2 = 3.6e-3;
    SW4 = 25.5e-3;
    
    % patch definitions
    patches.Iron = [2,3,4,6,8,12] ;
    patches.Air = [1,5,9,10];
    patches.Windings = [7,11];
    
    if nargin == 1
        param_names = fieldnames (parameters);
        for iParam  = 1:numel (param_names)
          eval ([param_names{iParam} '= parameters.(param_names{iParam});']);
        end
    end
    % calculated parameters
    GAP = (SD2 - RD1)/2;
    parameters.OUTD = SD1;
    alpha = pi/SLOTS;

    a1 = asin((SW1/2)/(SD2/2));
    a2 = asin((SW1/2)/(SD2/2+ST));
    a3 = asin((SW1/2)/(SD1/2-SW4));
    a4 = 2*pi/SLOTS - 2*a1;
    a5 = 2*pi/SLOTS - 2*a2;
    a6 = 2*pi/SLOTS - 2*a3;
    a7 = 2*asin(SW2/SD2);
    a8 = 2*asin((SW2/2)/(SD2/2+ST));

    srf(1) = nrbruled(nrbcirc(SD2/2-GAP/2,[0,0],0,a1),nrbcirc(SD2/2,[0,0],0,a1));
    srf(2) = nrbruled(nrbcirc(SD2/2,[0,0],0,a1),nrbcirc(SD2/2+ST,[0,0],0,a2));
    srf(3) = nrbruled(nrbcirc(SD2/2+ST,[0,0],0,a2),nrbcirc(SD1/2-SW4,[0,0],0,a3));
    srf(4) = nrbruled(nrbcirc(SD1/2-SW4,[0,0],0,a3),nrbcirc(SD1/2,[0,0],0,a3));

    srf(5) = nrbruled(nrbcirc(SD2/2-GAP/2,[0,0],a1,a1+(a4-a7)/2),nrbcirc(SD2/2,[0,0],a1,a1+(a4-a7)/2));
    srf(6) = nrbruled(nrbcirc(SD2/2,[0,0],a1,a1+(a4-a7)/2),nrbcirc(SD2/2+ST,[0,0],a2,a2+(a5-a8)/2));
    srf(7) = nrbruled(nrbcirc(SD2/2+ST,[0,0],a2,a2+(a5-a8)/2),nrbcirc(SD1/2-SW4,[0,0],a3,a3+(a6-a8)/2));
    srf(8) = nrbruled(nrbcirc(SD1/2-SW4,[0,0],a3,a3+(a6-a8)/2),nrbcirc(SD1/2,[0,0],a3,a3+(a6-a8)/2));

    srf(9) = nrbruled(nrbcirc(SD2/2-GAP/2,[0,0],a1+(a4-a7)/2,a1+(a4-a7)/2+a7/2),nrbcirc(SD2/2,[0,0],a1+(a4-a7)/2,a1+(a4-a7)/2+a7/2));
    srf(10) = nrbruled(nrbcirc(SD2/2,[0,0],a1+(a4-a7)/2,a1+(a4-a7)/2+a7/2),nrbcirc(SD2/2+ST,[0,0],a2+(a5-a8)/2,a2+(a5-a8)/2+a8/2));
    srf(11) = nrbruled(nrbcirc(SD2/2+ST,[0,0],a2+(a5-a8)/2,a2+(a5-a8)/2+a8/2),nrbcirc(SD1/2-SW4,[0,0],a3+(a6-a8)/2,a3+(a6-a8)/2+a8/2));
    srf(12) = nrbruled(nrbcirc(SD1/2-SW4,[0,0],a3+(a6-a8)/2,a3+(a6-a8)/2+a8/2),nrbcirc(SD1/2,[0,0],a3+(a6-a8)/2,a3+(a6-a8)/2+a8/2));
    % mirror first half winding
    nsrf = numel(srf);
    mirror = diag(ones(4,1));
    mirror(2,2) = -1;
    for i = 1:nsrf
        srf(nsrf+i) = nrbtform(srf(i), mirror*vecrotz(-(2*pi/SLOTS)));
    end
    patches.Iron = [patches.Iron, patches.Iron+nsrf];
    patches.Air = [patches.Air, patches.Air+nsrf];
    patches.Windings = [patches.Windings, patches.Windings+nsrf];
    % mirror rest
    nsrf = numel(srf);
    IronInit = patches.Iron;
    AirInit = patches.Air;
    WindingsInit = patches.Windings;
    for I = 1:5
        for i = 1:nsrf
            srf(I*nsrf+i) = nrbtform(srf(i), mirror*vecrotz(-(2*pi/SLOTS*(I+1))));
        end
        patches.Iron = [patches.Iron, IronInit+nsrf*I];
        patches.Air = [patches.Air, AirInit+nsrf*I];
        patches.Windings = [patches.Windings, WindingsInit+nsrf*I];
    end

    if (draw_geometry)
        figure(43)
        clf
        hold on
        for i = 1:numel(srf)
            if ismember(i, patches.Iron)
                nrbplotcol(srf(i), [10,10], 'color', [0.3,0.3,0.3]);
            elseif ismember(i, patches.Air)
                nrbplotcol(srf(i), [10,10], 'color', [0.1,0.1,1]);
            elseif ismember(i, patches.Windings)
                nrbplotcol(srf(i), [10,10], 'color', [0.9,0.1,1]);
            else
                disp("Stator patch is not material defined")
                disp(i)
                nrbplotcol(srf(i), [10,10], 'color', [0.1,0.1,0.1]);
            end
        end

        view(2)
        axis equal
    end

end