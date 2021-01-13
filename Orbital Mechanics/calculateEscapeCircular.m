%% calculateEscapeCircular.m
%Returns a variety of information about possible trajectories using a given
%planet as a stepping stone to escape the solar system
%inputs:
%planetid, the number of the planet you wish to use in the transfer (1-8,
%excludinig 3)
%
%outputs:
%struct DV, contains:
%   ... (TODO: fill in output)


function [DV] = calculateEscapeCircular(id)
    %% Circular approximation bookkeeping
    Re = 6378;
    radP = [2439; %mercury
            6051.85; %venus
            6378; %earth
            3396; %mars
            69911; %jupiter
            9*Re; %saturn
            3.93*Re; %uranus
            3.87*Re]; %Neptune

    flybyPeriapseMin = 1.025*radP;

    ue = 398600;
    us = 1.327e11;
    planetMu = [.055*ue; %mercury
                .815*ue; %venus
                ue; %earth
                .107*ue; %mars
                318*ue; %Jupiter
                95*ue; %saturn
                14.5*ue;%uranus
                17.2*ue]; %neptune

    rc = [  5.79e7; %mercury
            1.08e8; %venus
            1.496e8; %earth
            2.28e8; %mars
            7.78e8; %jupiter
            1.43e9; %saturn
            2.87e9; %uranus
            4.5e9]; %neptune

    %calculate approximate vinf coming from earth transfer
    deltaV1 = sqrt(us./rc(3)).*(sqrt(2.*rc./(rc(3) + rc))-1);
    vinf = sqrt(us./rc).*(1-sqrt(2.*rc(3)./(rc(3) + rc)));

    
    rPH = [rc(id);0;0];
    vPH= [0;sqrt(us/rc(id));0];
    vSH = vPH - [0;vinf(id);0];
    d = Flyby(us,planetMu(id),rPH,vPH,vSH,flybyPeriapseMin(id),'dark',flybyPeriapseMin(id),'Max Turn Hohmann',1);
    d.delta 
    DV.initialTransferDV = deltaV1(id);
    DV.escapeAtPlanetDark = sqrt(2*us/norm(rPH));
    DV.dvEscapeAtPlanetDark = DV.escapeAtPlanetDark - norm(vSH + [d.dVH;0]);
    DV.dvTotalAtPlanetDark = abs(deltaV1(id)) + max(0,DV.dvEscapeAtPlanetDark);
    DV.rPSHDark = abs(d.aF*(1-d.eF));
    DV.rASHDark = abs(d.aF*(1+d.eF));
    DV.escapePeriapseDark = sqrt(2*us/norm(DV.rPSHDark));
    DV.vPeriapseDark = us/d.hF*(1+d.eF);
    DV.vApoapsisDark = us/d.hF*(1-d.eF);
    DV.dvEscapePeriapseDark = DV.escapePeriapseDark - DV.vPeriapseDark;
    DV.dvTotalAtPeriapseDark = abs(deltaV1(id)) + max(0,DV.dvEscapePeriapseDark);
    DV.eF = d.eF;
    %note: repeating for light side doesn't matter cause symmetry (just
    %changes where on new orbit you end up and its orientation
    %repeat for light side
%     d = Flyby(us,planetMu(id),rPH,vPH,vSH,flybyPeriapseMin(id),'light',flybyPeriapseMin(id),'Max Turn Hohmann',1);
%     DV.initialTransferDV = deltaV1(id);
%     DV.escapeAtPlanetLight = sqrt(2*us/norm(rPH));
%     DV.dvEscapeAtPlanetLight = DV.escapeAtPlanetLight - norm(vSH + [d.dVH;0]);
%     DV.dvTotalAtPlanetLight = abs(deltaV1(id)) + max(0,DV.dvEscapeAtPlanetLight);
%     DV.rSHLight = abs(d.aF*(1-d.eF));
%     DV.escapePeriapseLight = sqrt(2*us/norm(DV.rSHLight));
%     DV.vPeriapseLight = us/d.hF*(1+d.eF);
%     DV.vApoapsisLight = us/d.hF*(1-d.eF);
%     DV.dvEscapePeriapseLight = DV.escapePeriapseLight - DV.vPeriapseLight;
%     DV.dvTotalAtPeriapseLight = abs(deltaV1(id)) + max(0,DV.dvEscapePeriapseLight)
    DV.d = d;
    DV
end