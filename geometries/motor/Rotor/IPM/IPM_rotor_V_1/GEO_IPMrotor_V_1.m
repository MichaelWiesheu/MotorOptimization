function [srf, patches, parameters] = GEO_IPMrotor_V_1 (parameters)
    draw_geometry = true;
    % drive standard parameters
    %GAP = 1e-3;
    HEIGHT = 35e-3;
    POLES = 4;
    SLOTS = 24;
    % rotor standard parameters
    DMAG = 30e-3;
    LSLIT1 = 6.4e-3;
    LSLIT2 = 4.3e-3;
    DSLIT5 = 1.0e-3;
    DSLIT6 = 2.0e-3;
    MA = 150;
    MT1 = 4.0e-3;
    MW1 = 22.0e-3;
    RA1 = 144.2;
    RA2 = 166.0;
    RD1 = 100e-3;
    RD2 = 29.1e-3;
    RF = 0.2e-3;
    RS = 1.0e-3;
    RW2 = 1.0e-3;
    RW3 = 1.0e-3;
    RW4 = 1.0e-3;
    RW5 = 1.0e-3;
    SD2 = 102e-3;
    WMAG = 4.0e-3;
    
    REFINEMENT = [1, 0, 3, 0, 1];

    % patch definitions
    patches.Iron = [2, 3, 7, 8, 12:14] ;
    patches.Magnets = [9];
    patches.Air = [1, 4, 5, 6, 10, 11];
    patches.Windings = [];
    % overwrite standard parameters by input data
    if nargin == 1
        param_names = fieldnames (parameters);
        for iParam  = 1:numel (param_names)
          eval ([param_names{iParam} '= parameters.(param_names{iParam});']);
        end
    end
    % calculated parameters
    GAP = (SD2 - RD1)/2;
    beta = pi/POLES;
    motorAngle = 2*beta;
    gamma = atan(WMAG/DMAG);
    delta = beta-gamma;
    epsilon = deg2rad(MA/2)-beta;
    phi1 = deg2rad(rad2deg(epsilon) + RA1 - 180);
    phi2 = deg2rad(rad2deg(epsilon) + RA2 - 180);

    % point definitions
    point{1} = [cos(beta), sin(beta)]*RD2/2;
    % magnet
    point{2} = [cos(delta), sin(delta)]*(DMAG^2+WMAG^2)^0.5;
    point{3} = point{2} + [cos(epsilon), -sin(epsilon)]*MW1;
    point{4} = point{2} + [-sin(epsilon), -cos(epsilon)]*MT1;
    point{5} = point{3} + [-sin(epsilon), -cos(epsilon)]*MT1;
    % rw points
    point{6} = point{2} + [-cos(epsilon), sin(epsilon)]*RW4; 
    point{7} = point{3} + [cos(epsilon), -sin(epsilon)]*RW2; 
    point{8} = point{4} + [-cos(epsilon), sin(epsilon)]*RW5; 
    point{9} = point{5} + [cos(epsilon), -sin(epsilon)]*RW3; 

    point{10} = [(DMAG-DSLIT5)*cos(beta), (DMAG-DSLIT5)*sin(beta)];
    point{11} = [(DMAG-DSLIT6)*cos(beta), (DMAG-DSLIT6)*sin(beta)];
    
    point{12} = point{10} + [sin(beta), -cos(beta)]*RS/2;
    point{13} = point{11} + [sin(beta), -cos(beta)]*RS/2;

    point{14} = point{9} + [cos(phi1), -sin(phi1)]*LSLIT1;
    point{15} = point{7} + [cos(phi2), -sin(phi2)]*LSLIT2;

    point{16} = [RD2/2, 0];
    point{17} = [point{5}(1), 0];
    point{18} = [RD1/2, 0];
    
    point{19} = vecnormalize(point{15})*RD1/2;
    point{20} = vecnormalize(point{3})*RD1/2;
    point{21} = vecnormalize(point{2})*RD1/2;
    point{22} = vecnormalize(point{13})*RD1/2;

    point{23} = point{10} - [sin(beta), -cos(beta)]*RS/2;
    point{24} = point{11} - [sin(beta), -cos(beta)]*RS/2;
    point{25} = vecnormalize(point{13})*RD2/2;

    % hard coded / vs projected
    a1 = motorAngle/9; 
    a2 = 2*motorAngle/9;  
    a3 = 3*motorAngle/9; 
    a4 = beta-atan(RS/RD1);
%     a1 = atan2(point{19}(2), point{19}(1));
%     a2 = atan2(point{20}(2), point{20}(1));
%     a3 = atan2(point{21}(2), point{21}(1));
%     a4 = atan2(point{22}(2), point{22}(1));
    
    % circs defining the refinement of the motor
    %circ1 = nrbcirc(RD1/2, [0, 0], 0, a1);
    circ1 = nrbcirc(RD1/2, [0, 0], 0, a1);
    circ2 = nrbcirc(RD1/2, [0, 0], a1, a2);
    circ3 = nrbcirc(RD1/2, [0, 0], a2, a3);
    circ4 = nrbcirc(RD1/2, [0, 0], a3, a4);
    circ5 = nrbcirc(RD1/2, [0, 0], a4, beta+beta-a4);
    % refine circle as wanted...
    circ1 = nrbrefine1D(circ1, REFINEMENT(1));
    circ2 = nrbrefine1D(circ2, REFINEMENT(2));
    circ3 = nrbrefine1D(circ3, REFINEMENT(3));
    circ4 = nrbrefine1D(circ4, REFINEMENT(4));
    circ5 = nrbrefine1D(circ5, REFINEMENT(5));

    % Hard coded values for air gap, which needs to remain the same
    A1 = motorAngle/9; 
    A2 = 2*motorAngle/9;
    A3 = 3*motorAngle/9; 
    A4 = 4*motorAngle/9;
    srf(1) = nrbruled(circ1, nrbcirc(RD1/2+GAP/2, [0, 0], 0, A1));
 
    % air right slit
    [RS_Top1, RS_1fix1, RS_Top2] = line2radius(nrbline(point{3}, point{7}), nrbline(point{7}, point{15}), RF);
    [RS_1fix2, RS_2fix1, RS_Right1] = line2radius(nrbline(point{7}, point{15}), nrbline(point{15}, point{14}), RF);
    [RS_2fix2, RS_3fix1, RS_Right3] = line2radius(nrbline(point{15}, point{14}), nrbline(point{14}, point{9}), RF);
    [RS_3fix2, RS_Bot3, RS_Bot2] = line2radius(nrbline(point{14}, point{9}), nrbline(point{9}, point{5}), RF);
    % If inner angle of right curves is > 90 deg, manually increase, such that
    % geometry can be compared to ones with angles of < 90 deg
    if numel(RS_Top2.knots) == 6
        RS_Top2 = nrbkntins(RS_Top2, [0.5, 0.5]);
    end
    if numel(RS_Right1.knots) == 6
        RS_Right1 = nrbkntins(RS_Right1, [0.5, 0.5]);
    end
    if numel(RS_Right3.knots) == 6
        RS_Right3 = nrbkntins(RS_Right3, [0.5, 0.5]);
    end
    if numel(RS_Bot2.knots) == 6
        RS_Bot2 = nrbkntins(RS_Bot2, [0.5, 0.5]);
    end

%     RS_Right1 = nrbdegelev(RS_Right1, 1);
%     RS_Right3 = nrbdegelev(RS_Right3, 1);
    RS_Top3 = nrbline(RS_1fix1.coefs(1:2, 1), RS_1fix2.coefs(1:2, 2));
    temp = nrbglue(nrbdegelev(RS_Top1, 1), RS_Top2, 2, 1);
    RS_Top = nrbglue(temp, nrbdegelev(RS_Top3, 1), 2, 1);
    RS_Top = nrbglue(RS_Top, RS_Right1, 2, 1);
    RS_Top.knots = RS_Top.knots / max(RS_Top.knots);
    
    RS_Right = nrbline(RS_2fix1.coefs(1:2, 1), RS_2fix2.coefs(1:2, 2));

    RS_Bot1 = nrbline(RS_3fix1.coefs(1:2, 1), RS_3fix2.coefs(1:2, 2));
    temp = nrbglue(RS_Right3, nrbdegelev(RS_Bot1, 1), 2, 1);
    temp = nrbglue(temp, RS_Bot2, 2, 1);
    RS_Bot = nrbglue(temp, nrbdegelev(RS_Bot3, 1), 2, 1);
    RS_Bot.knots = RS_Bot.knots / max(RS_Bot.knots);

    % right iron
    srf(2) = nrbruled(nrbreverse(RS_Right), circ1);
    % air
    srf(3) = nrbruled(nrbreverse(RS_Top), circ2);
    % air slit right
    srf(4) = nrbcoons(nrbreverse(RS_Bot), nrbreverse(nrbextract(srf(3), 3)), nrbline(point{5}, point{3}), nrbextract(srf(2), 3));
    % air
    srf(5) = nrbruled(nrbextract(srf(3), 4), nrbcirc(RD1/2 + GAP/2, [0, 0], A1, A2));
    srf(6) = nrbruled(circ3, nrbcirc(RD1/2 + GAP/2, [0, 0], A2, A3));
    % iron
    srf(7) = nrbruled(circ3, nrbline(point{3}, point{2}));

    % Air slit left
    [LS_Top1, LS_1fix1, LS_Top2] = line2radius(nrbline(point{2}, point{6}), nrbline(point{6}, point{12}), RF);
    [LS_1fix2, LS_2fix1, LS_Left1] = line2radius(nrbline(point{6}, point{12}), nrbline(point{12}, point{13}), RF);
    [LS_2fix2, LS_3fix1, LS_Left3] = line2radius(nrbline(point{12}, point{13}), nrbline(point{13}, point{8}), RF);
    [LS_3fix2, LS_Bot3, LS_Bot2] = line2radius(nrbline(point{13}, point{8}), nrbline(point{8}, point{4}), RF);
    if numel(LS_Top2.knots) == 6
        LS_Top2 = nrbkntins(LS_Top2, [0.5, 0.5]);
    end
    if numel(LS_Left1.knots) == 6
        LS_Left1 = nrbkntins(LS_Left1, [0.5, 0.5]);
    end
    if numel(LS_Left3.knots) == 6
        LS_Left3 = nrbkntins(LS_Left3, [0.5, 0.5]);
    end
    if numel(LS_Bot2.knots) == 6
        LS_Bot2 = nrbkntins(LS_Bot2, [0.5, 0.5]);
    end

    LS_Top3 = nrbline(LS_1fix1.coefs(1:2, 1), LS_1fix2.coefs(1:2, 2));
    temp = nrbglue(nrbdegelev(LS_Top1, 1), LS_Top2, 2, 1);
    LS_Top = nrbglue(temp, nrbdegelev(LS_Top3, 1), 2, 1);
    LS_Top = nrbglue(LS_Top, LS_Left1, 2, 1);
    LS_Top.knots = LS_Top.knots / max(LS_Top.knots);
    
    LS_Left = nrbline(LS_2fix1.coefs(1:2, 1), LS_2fix2.coefs(1:2, 2));

    LS_Bot1 = nrbline(LS_3fix1.coefs(1:2, 1), LS_3fix2.coefs(1:2, 2));
    temp = nrbglue(LS_Left3, nrbdegelev(LS_Bot1, 1), 2, 1);
    temp = nrbglue(temp, LS_Bot2, 2, 1);
    LS_Bot = nrbglue(temp, nrbdegelev(LS_Bot3, 1), 2, 1);
    LS_Bot.knots = LS_Bot.knots / max(LS_Bot.knots);
    % iron
    srf(8) = nrbruled(circ4, LS_Top);
    % magnet
    srf(9) = nrbcoons(nrbline(point{4}, point{5}), nrbreverse(nrbextract(srf(7), 4)), nrbline(point{4}, point{2}), nrbextract(srf(4), 1));
    % air slit left
    srf(10) = nrbcoons(nrbreverse(LS_Bot), nrbextract(srf(8), 4), nrbextract(srf(9), 1), nrbreverse(LS_Left));
    % air 
    srf(11) = nrbruled(nrbextract(srf(8), 3), nrbcirc(RD1/2 + GAP/2, [0, 0], A3, A4));
    % iron
    srf(12) = nrbruled(nrbextract(srf(9), 3), nrbline(point{16}, point{17}));
    srf(13) = nrbruled(nrbextract(srf(4), 3), nrbline(point{17}, point{18}));%nrbcirc(RD1/2, [0, 0], 0, a1));
    srf(14) = nrbruled(nrbcirc(RD2/2, [0, 0], 0, a4), nrbextract(srf(10), 3));
    % in between patches
    p1left = LS_Left.coefs(1:2, 1) + [-sin(beta); cos(beta)]*RS;
    p2left = LS_Left.coefs(1:2, 2) + [-sin(beta); cos(beta)]*RS;
    srf(15) = nrbruled(circ5, nrbline(LS_Left.coefs(1:2, 1), p1left));
    srf(16) = nrbcoons(nrbextract(srf(10), 2), nrbline(p2left, p1left), nrbline(LS_Left.coefs(1:2, 2), p2left), nrbextract(srf(15), 4));
    srf(17) = nrbruled(nrbextract(srf(16), 1), nrbcirc(RD2/2, [0, 0], a4, beta+(beta-a4)));
    srf(18) = nrbruled(circ5, nrbcirc(RD1/2 + GAP/2, [0, 0], A4, beta+(beta-A4)));
    % mirror patches
    nsrf = numel(srf);
    mirror = diag(ones(4, 1));
    mirror(2, 2) = -1;
    for i = 1:14
        srf(nsrf+i) = nrbtform(srf(i), mirror*vecrotz(-2*beta));
    end
    patches.Iron = [patches.Iron, patches.Iron+nsrf, 15:17];
    patches.Air = [patches.Air, patches.Air+nsrf, 18];
    patches.Magnets = [patches.Magnets, patches.Magnets+nsrf];


    if (draw_geometry)
        figure(42)
        clf
        %line([0, 0.1], [0, 0.1])
        axis equal
        hold on
        for i = 1:numel(srf)
            if ismember(i, patches.Iron)
                nrbplotcol(srf(i), [10, 10], 'color', [0.3, 0.3, 0.3]);
            elseif ismember(i, patches.Magnets)
                nrbplotcol(srf(i), [10, 10], 'color', 'green');
            elseif ismember(i, patches.Air)
                nrbplotcol(srf(i), [10, 10], 'color', [0.1, 0.1, 1]);
            elseif ismember(i, patches.Windings)
                nrbplotcol(srf(i), [10, 10], 'color', [0.9, 0.1, 1]);
            else
                disp("Rotor patch is not material defined")
                disp(i)
                nrbplotcol(srf(i), [10, 10], 'color', [0.1, 0.1, 0.1]);
            end
        end

        for i = 1:numel(point)
            scatter(point{i}(1), point{i}(2), "black", "filled");
            text(point{i}(1)+1e-4, point{i}(2), string(i))
        end
        view(2)
        axis equal
    end

end

function [l1, l2, arc] = line2radius(line1, line2, radius)
    assert(size(line1.coefs, 2)==2, "Line1 needs to be a line");
    assert(size(line2.coefs, 2)==2, "Line2 needs to be a line");
    assert(all(line1.coefs(1:2, 2)==line2.coefs(1:2, 1)), "Endpoint of line1 needs to be starting point of line2");

    p1 = line1.coefs(1:2, 1);
    p2 = line1.coefs(1:2, 2);
    p3 = line2.coefs(1:2, 2);

    sp1 = vecnormalize(line1.coefs(1:2, end)-line1.coefs(1:2, 1));
    sp2 = vecnormalize(line2.coefs(1:2, end)-line2.coefs(1:2, 1));
    
    alpha = acos(sum(sp1.*sp2));
    % If angle is too small, replace radius by line
    if alpha < 1e-4
        warning("JMAG_rotor: Small radius angle! Radius is removed!");
        frac = 0.001; % share of "arc" which is replaced by line
        l1 = nrbline(p1, frac*p1 + (1-frac)*p2);
        arc = nrbline(frac*p1 + (1-frac)*p2, frac*p3 + (1-frac)*p2);
        arc = nrbdegelev(arc, 1);
        l2 = nrbline(frac*p3 + (1-frac)*p2, p3);
        return
    end
    d = radius*tan(alpha/2);
    P1 = p2-d*sp1;
    l1 = nrbline(p1, P1);
    P2 = p2+d*sp2;
    l2 = nrbline(P2, p3);
    rvec = vecnormalize(0.5*(sp1-sp2));
    Pr = p2-(radius^2+d^2)^0.5*rvec;
    
    angle1 = atan2(P1(2)-Pr(2), P1(1)-Pr(1));
    angle2 = atan2(P2(2)-Pr(2), P2(1)-Pr(1));
    
    angle1 = angle1 + 2*pi*(angle1<0);
    angle2 = angle2 + 2*pi*(angle2<0);
    if angle1 > angle2
        angle1 = angle1 - 2*pi;
    end
    if angle2 > angle1 + pi
        temp = angle1;
        angle1 = angle2;
        angle2 = temp;
        arc = nrbreverse(nrbcirc(radius, Pr, angle1, angle2));
    else
        arc = nrbcirc(radius, Pr, angle1, angle2);
    end
end

function crv = nrbrefine1D (crv, nctrlptsins) % switch to kntrefine?
    vals = linspace(0, 1, nctrlptsins+2);
    vals(nctrlptsins+2) = [];
    vals(1) = [];
    crv = nrbkntins(crv, vals);
end