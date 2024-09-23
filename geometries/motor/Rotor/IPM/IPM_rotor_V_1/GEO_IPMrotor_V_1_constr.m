function c = GEO_IPMrotor_V_1_constr (parameters)
    draw_geometry = true;
    % drive standard parameters
    GAP = 1e-3;
    % GapType = 'free'  TBD usage in JMAG?
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

    SD1 = 204e-3;
    ST = 1.6e-3;
    SW1 = 6.0e-3;
    SW2 = 3.6e-3;
    SW4 = 25.5e-3;

    % overwrite standard parameters by input data
    if nargin == 1
        param_names = fieldnames (parameters);
        for iParam  = 1:numel (param_names)
          eval ([param_names{iParam} '= parameters.(param_names{iParam});']);
        end
    end
    % calculated parameters
    beta = pi/POLES;
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

    point{16} = [RD2/2,0];
    point{17} = [point{5}(1),0];
    point{18} = [RD1/2,0];
    
    point{19} = vecnormalize(point{15})*RD1/2;
    point{20} = vecnormalize(point{3})*RD1/2;
    point{21} = vecnormalize(point{2})*RD1/2;
    point{22} = vecnormalize(point{13})*RD1/2;

    point{23} = point{10} - [sin(beta), -cos(beta)]*RS/2;
    point{24} = point{11} - [sin(beta), -cos(beta)]*RS/2;
    point{25} = vecnormalize(point{13})*RD2/2;

    minIron = 1.5e-3;
    c(1) = vecmag(point{9}) - (RD1/2) + minIron;
    c(2) = vecmag(point{14}) - (RD1/2) + minIron;
    c(3) = -point{14}(2) + minIron/2;
    c(4) = -point{9}(2) + minIron/2;
    c(5) = vecmag(point{15}) - (RD1/2) + minIron;
    c(6) = -vecmag(point{15}-point{14}) + 2*RF;
    c(7) = -vecmag(point{12}-point{13}) + 2*RF;
    c(8) = RD2/2 + 2*minIron -vecmag(point{8}) ;
    c(9) = RD2/2 + 2*minIron -vecmag(point{13});

    % Stator:
    % a1 = asin((SW1/2)/(SD2/2));
    % a7 = 2*asin(SW2/SD2);
    % c(10) = -(pi/SLOTS-a1-0.5*a7);
    % c(11) = -(SD1/2 - SW4 - ST - SD2/2);

    c = c * 100;    % scale to create better convergence behavior


end

