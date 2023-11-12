classdef nozzleGeometry
    %obj = nozzleGeometry(expRatio_in, contRatio_in, Dt_in,...
    %rCon_in, thetaCon_in, rExp_in, thetaExp_in, Lc_in)
    %Creates a nozzleGeometry object with the following properties:
    %ExpRatio = Expansion ratio
    %Dc: Chamber diameter [m]
    %Dt: Throat diameter [m]
    %rCon: Chamber-to-Contraction curve radius [m]
    %thetaCon: Angle between chamber wall and contracting nozzle wall [deg]
    %rExp: Contraction-to-expansion curve radius [m]
    %thetaExp: Angle between throat axis and expanding nozzle wall [deg]
    %Lc: Chamber length [m]
    
    properties
    %Diameters and area ratios
        expRatio
        contRatio
        Dt
        
        %Contraction and Expansion Geometries
        rCon
        rExp
        thetaCon
        thetaExp
        
        %Nozzle section lengths
        Lc
        LCon
        LExp
    
        %% Calculated properties
        %Diameters
        De
        Dc
        D2
        D3
        D5
        
        %Lengths
        L2
        L3
        Lt
        L5
        L_total
        
        %Areas
        At
        Ae
        Ac
    end
    
    methods
        function obj = nozzleGeometry(expRatio_in, contRatio_in, Dt_in, rCon_in,...
                                      thetaCon_in, rExp_in, thetaExp_in,...
                                      Lc_in)
            %%NOTE:
            %%Geometries are collected for the nozzle as required by 
            %%USAF tech. report AEDC-TR-91-1
            
            %Define given properties using constructor inputs
            obj.expRatio = expRatio_in;
            obj.contRatio = contRatio_in;
            obj.Dt = Dt_in;
            
            obj.rCon = rCon_in;
            obj.thetaCon = thetaCon_in;
            obj.rExp = rExp_in;
            obj.thetaExp = thetaExp_in;
            
            obj.Lc = Lc_in;
            
            
            %Calculate dependent properties using given arguments
            obj.Dc = obj.Dt*sqrt(obj.contRatio);
            obj.De = obj.Dt*sqrt(obj.expRatio);
            obj.D2 = obj.Dc - 2*obj.rCon*(1-cosd(obj.thetaCon));
            obj.D3 = obj.Dt + 2*obj.rExp*(1-cosd(obj.thetaCon));
            obj.D5 = obj.Dt + 2*obj.rExp*(1-cosd(obj.thetaExp));
            obj.L2 = obj.Lc + obj.rCon*sind(obj.thetaCon);
            obj.L3 = obj.L2 + (obj.D2 - obj.D3)/(2*tand(obj.thetaCon));
            obj.Lt = obj.L3 + obj.rExp*sind(obj.thetaCon);
            obj.L5 = obj.Lt + obj.rExp*sind(obj.thetaExp);
            obj.L_total = obj.L5 + (obj.De-obj.D5)/(2*tand(obj.thetaExp));
            
            %Calculate throat area
            obj.At = 0.25*pi*obj.Dt.^2;
            obj.Ae = obj.At*obj.expRatio;
            obj.Ac = 0.25*pi*obj.Dc.^2;
        end
    end
end

