classdef injectorGeometry < handle
    %INJECTORGEOMETRY Injector with given orifice geometry
    %
    %Constructor format:
    %
    %   object = injectorGeometry(oxList, fuelList, BLCList)
    %
    %Element list format:
    %
    %   The list for each injection type is formatted as follows:
    %
    %   [Cd1, A1, count1; Cd2, A2, count2;...]
    %
    %   So for example, if oxList is equal to...
    %
    %   oxList == [0.7, 1E-6, 24; 0.8, 3E-6, 16]
    %
    %   ...then the injector object will have 24 oxidizer inlets with
    %   discharge coefficients of 0.7 and outlet areas of 1E-6 and 16
    %   oxidizer inlets with discharge coefficients of 0.8 and outlet areas
    %   of 3E-6.
    %
    %   If the list is left empty then the CdAs will be 0.
    
    properties
        %%Discharge coefficients for each injection type
        %Average discharge coefficient of all ox inlets
        oxCd
        %Average discharge coefficient of all fuel inlets
        fuelCd
        %Average discharge coefficeint of all boundary layer coolant inlets
        BLCCd
        
        %%Total inlet areas for each type across the injector
        %Injector total ox inlet area
        oxArea
        %Injector total fuel inlet area
        fuelArea
        %Injector total boundary layer coolant area
        BLCArea
        
    end
    
    methods
        function obj = injectorGeometry(oxList, fuelList, BLCList)
            %INJECTORGEOMETRY Create an injector with specified geometry.
            %
            %   obj = injectorGeometry(oxList, fuelList, BLCList)
            %
            %Element list format:
            %
            %   The list for each injection type is formatted as follows:
            %
            %   [Cd1, A1, count1; Cd2, A2, count2;...]
            %
            %   So for example, if oxList is equal to...
            %
            %   oxList == [0.7, 1E-6, 24; 0.8, 3E-6, 16]
            %
            %   ...then the injector object will have 24 oxidizer inlets with
            %   discharge coefficients of 0.7 and outlet areas of 1E-6 and 16
            %   oxidizer inlets with discharge coefficients of 0.8 and outlet areas
            %   of 3E-6.
            %
            %   If the list is left empty then the CdAs will be 0.
            
            %Initialize values at zero
            obj.oxCd = 0;
            obj.oxArea = 0;
            
            if nargin > 0
                %Number of oxidizer inlet types, number of rows in oxList
                numOxTypes = size(oxList, 1);
                for i = 1:numOxTypes
                    Cdi = oxList(i, 1);
                    Ai = oxList(i, 2);
                    counti = oxList(i, 3);

                    %Take area-weighted average of Cd across injector for
                    %given injection type
                    obj.oxCd = (obj.oxCd*obj.oxArea + Cdi*Ai*counti)/(obj.oxArea + Ai*counti);
                    obj.oxArea = obj.oxArea + Ai*counti;

                    %Repeat until all rows/element types have been lumped in
                end
            end
                
            %Repeat for other injection types
            obj.fuelCd = 0;
            obj.fuelArea = 0;
            
            if nargin > 1
                numFuelTypes = size(fuelList, 1);
                for i = 1:numFuelTypes
                    Cdi = fuelList(i, 1);
                    Ai = fuelList(i, 2);
                    counti = fuelList(i, 3);

                    obj.fuelCd = (obj.fuelCd*obj.fuelArea + Cdi*Ai*counti)/(obj.fuelArea + Ai*counti);
                    obj.fuelArea = obj.fuelArea + Ai*counti;
                end
            end
            
            obj.BLCCd = 0;
            obj.BLCArea = 0;
            
            if nargin > 2
                numBLCTypes = size(BLCList, 1);
                for i = 1:numBLCTypes

                    Cdi = BLCList(i, 1);
                    Ai = BLCList(i, 2);
                    counti = BLCList(i, 3);

                    obj.BLCCd = (obj.BLCCd*obj.BLCArea + Cdi*Ai*counti)/(obj.BLCArea + Ai*counti);
                    obj.BLCArea = obj.BLCArea + Ai*counti;
                end
            end 
            
        end
            
        %%Add Element Methods
        %These methods can add additional orifices. Likely won't see much
        %use if it's made right the first time, but the feature itself is
        %nice to have.
        
        function [] = addOx(obj, Cd, area, count)
            %   object = object.addOx(Cd, Area, count);
            %
            %   Adds (count) number of oxidizer injection elements with
            %   given CdA 
            
                %Get the old properties...
                oldCd = obj.oxCd;
                oldA  = obj.oxArea;
                
                %...average them with the new properties to get the average
                %across all ox elements
                obj.oxCd = (oldA*oldCd + count*Cd*area)/(oldA + count*area);
                obj.oxArea = oldA + count*area;
        end
        
        function addFuel(obj, Cd, area, count)
            %   object = object.addFuel(Cd, Area, count);
            %
            %   Adds (count) number of fuel injection elements with
            %   given CdA 
            
                %Get the old properties...
                oldCd = obj.fuelCd;
                oldA  = obj.fuelArea;
                
                %...average them with the new properties to get the average
                %across all ox elements
                obj.fuelCd = (oldA*oldCd + count*Cd*area)/(oldA + count*area);
                obj.fuelArea = oldA + count*area;
                end
        
        function addBLC(obj, Cd, area, count)
            %   object = object.addOx(Cd, Area, count);
            %
            %   adds (count) number of oxidizer injection elements with
            %   given CdA 
            
                %Get the old properties...
                oldCd = obj.BLCCd;
                oldA  = obj.BLCArea;
                
                %...average them with the new properties to get the average
                %across all ox elements
                obj.BLCCd = (oldA*oldCd + count*Cd*area)/(oldA + count*area);
                obj.BLCArea = oldA + count*area;
        end
    end
end

