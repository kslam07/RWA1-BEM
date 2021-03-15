classdef BEMSolverNREL
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    % user-defined
    properties
        rho(1,1) double {mustBeNonnegative} = 1.225;      % density kg/m2
        uInf (1,1) double {mustBeNonnegative} = 10.0;     % m/s
        TSR (1,1) double {mustBeNonnegative} = 6.0;       % tip speed ratio
        rRotor (1,1) double {mustBeNonnegative} = 50.0;   % radius of rotor
        rTipRatio (1,1) double {mustBeNonnegative} = 1.0  % r_tip/R_rotor
        rRootRatio (1,1) double {mustBeNonnegative} = 0.2;% r_root/R_rotor
        nAnnulus (1,1) double {mustBeNonnegative} = 10;   % blade segments
        spacing (1,:) char {mustBeMember(spacing, ["0", "1", "cosine"])}...
                = '0';                                    % segment spacing
        nBlades (1,1) double {mustBeNonnegative} = 3;     % blades /rotors
        bladePitch (1,1) double = -2;                     % deg
        atol (1,1) double = 1e-6;
        nIter (1,1) uint16 = 100; 
        yawAngle (1,1) double = 0;
        nPsi (1,1) uint16 = 10;
    end
    
    % computed properties
    properties%(Access = )
        % rows: annulus; column: elements in segments (varying azimuth)
        locSegment (1, :) double;   % radial length of each blade seg.
        rR (:, 1) double;           % centre loc. segment in radial dir.
        psiSegment (1, :) double;   % azimuth angle of element in annulus
        alpha (:, :) double;        % local angle-of-attack
        phi (:, :) double;          % inflow angle
        aprime (:, :) double;       % tangential flow factor
        a (:, :) double;            % axial flow factor
        Az (:, :) double;           % azimuthal force
        Ax (:, :) double;           % axial force
        Cl (:, :) double;           % sectional lift coefficient
        Cd (:, :) double;           % sectional drag coefficient
        CT (:, :) double;           % thrust coefficients
        CN (:, :) double;           % normal force coefficients
        Cq (:, :) double;           % torque coefficient
        thrustIter (:, :, :) double;   % thrust iterations
        fTot (:, :) double;         % Prandtl correction
        gamma (:, :) double;        % Circulation
        areaAnnulus (:, 1) double   % area of each annuli
        
    end
    
    properties%(Access = protected)
        airfoilFilename = "polar_DU95W180.xlsx";
        fCL
        fCD
        amax
        amin
        Omega
    end
    
    methods
        function obj = init(obj)
            % initialize the class object
            obj         = createPolarSplines(obj);
            obj         = computeSegments(obj);
            obj         = computePsiArr(obj);
            obj.Omega   = obj.TSR * obj.uInf / obj.rRotor;
            obj.yawAngle = deg2rad(obj.yawAngle);
            
            % allocate memory for rest of matrices holding solutions
            obj.alpha   = zeros(obj.nAnnulus, obj.nPsi);
            obj.phi     = zeros(obj.nAnnulus, obj.nPsi);
            obj.aprime  = zeros(obj.nAnnulus, obj.nPsi);
            obj.a       = zeros(obj.nAnnulus, obj.nPsi);
            obj.Az      = zeros(obj.nAnnulus, obj.nPsi);
            obj.Ax      = zeros(obj.nAnnulus, obj.nPsi);
            obj.CT      = zeros(obj.nAnnulus, obj.nPsi);
            obj.CN      = zeros(obj.nAnnulus, obj.nPsi);
            obj.Cq      = zeros(obj.nAnnulus, obj.nPsi);
            obj.rR      = zeros(1, obj.nAnnulus);
            obj.thrustIter = zeros(obj.nAnnulus, obj.nPsi, obj.nIter);
        end
        
        function obj = computePsiArr(obj)
            psiArr = linspace(0, 2*pi, obj.nPsi);
            obj.psiSegment = psiArr;
        end
        
        function obj = computeSegments(obj)
            % compute segment lengths based on user-input spacing arg.
            % spacing == '0':
            % The segments are divided linearly
            % spacing == 'cosine':
            % The segments are cosine space divided
            
            % Due to rotor hub, the blade starts at:
            rStart = obj.rRootRatio;  % r/R
            rEnd = obj.rTipRatio;  % r/R
            
            if obj.spacing == '0'    % linear spacing
                segmentBounds = linspace(rStart, rEnd, obj.nAnnulus+1);
%                 lengthSegments = diff(segmentBounds);
                
            elseif obj.spacing == "cosine"
                segmentBounds = cosspace(rStart, rEnd, obj.nAnnulus+1);
%                 lengthSegments = diff(segmentBounds);
                
            else
                disp("unrecognized spacing, set to default")
                segmentBounds = linspace(rStart, rEnd, obj.nAnnulus+1);
%                 lengthSegments = diff(segmentBounds);
                
            end
            obj.locSegment = segmentBounds;
        end
        
        function obj = createPolarSplines(obj)
            % Create fitted functions from the airfoil polar data
            % Create two spline functions to fit the CL and CD data from
            % the polar data given to the class.
            data = table2array(readtable(obj.airfoilFilename));
            alphaRad = deg2rad(data(:,1));
            obj.fCL = fit(alphaRad, data(:, 2), "cubicinterp");
            obj.fCD = fit(alphaRad, data(:, 3), "cubicinterp");
            obj.amax = max(alphaRad);
            obj.amin = min(alphaRad);
        end

        function aSkew = skewWakeCorr(obj, aSegment)
            chi = obj.yawAngle.*(1+0.6.*aSegment);
            aSkew = aSegment.*(1+15*pi/64.*obj.rR.*tan(chi/2).*sin(obj.psiSegment));
        end
        
        function fTot = calcPrandtlCorr(obj, phi)
            rTip = obj.rRotor;
            rRoot = obj.rRootRatio * obj.rRotor;
            r = obj.rR .* obj.rRotor;
            temp1 = obj.nBlades/2 * (rTip - r) ./ (r.*abs(sin(phi)));
            fTip = 2/pi * acos(exp(-temp1));
            temp1 = obj.nBlades/2 * (r-rRoot) ./ (rRoot * abs(sin(phi)));
            fRoot = 2/pi * acos(exp(-temp1));
            fTot = fTip .* fRoot; 
        end

        function [cAx, cAz, clSegment, cdSegment, alphaSegment] ...
                = computeLoadsSegment(obj, uRotor, uTan, twistAngle)   
            
            theta = deg2rad(obj.bladePitch);                % blade pitch
            phiSegment = atan2(uRotor, uTan);               % flow angle
            alphaSegment = phiSegment - twistAngle - theta; % AoA
            
            % compute Cl and Cd for all elements in annulus
            clSegment = reshape(obj.fCL(alphaSegment), ...
                [obj.nAnnulus, obj.nPsi]);
            cdSegment = reshape(obj.fCD(alphaSegment), ...
                [obj.nAnnulus, obj.nPsi]);
            
            % take CL and CD from highest alpha in data if alpha is higher
            clSegment(alphaSegment > obj.amax) = obj.fCL(obj.amax);
            cdSegment(alphaSegment > obj.amax) = obj.fCD(obj.amax);
            
            % take CL and CD from lowest alpha in data if alpha is lower
            clSegment(alphaSegment < obj.amin) = obj.fCL(obj.amin);
            cdSegment(alphaSegment < obj.amin) = obj.fCD(obj.amin);
            
            % compute the normal and tangential coefficients
            cAz = clSegment.*sin(phiSegment) - cdSegment.*cos(phiSegment);
            cAx = clSegment.*cos(phiSegment) + cdSegment.*sin(phiSegment);
        end
        
        function obj = solveStreamtube(obj)
            % solve balance of momentum between blade element load and
            % loading in the streamtube
            % rR1      : location of inner boundary of blade segment 
            % rR2      : location of outer boundary of blade segment
            % rRoot    : root of blade/rotor location in frac. of radius
            % rTip     : tip of blade/rotor location in fraction of radius
            % rRotor   : total radius of the rotor
            % uInf     : freestream velocity
            % Omega    : rotational velocity
            % nBlades  : number of rotor blades
            
            % center of blade segment
            obj.rR = (obj.locSegment(1:end-1)+ obj.locSegment(2:end))/2;
            % twist of blade seg.
            twistAngle = obj.computeTwists(obj.rR);
            % chord length
            chord = obj.computeChordLength(obj.rR); 
            % area of blade elements with shape (element x annulus)
            areaSegment = (pi * ((obj.locSegment(2:end)*obj.rRotor).^2 ...
                - (obj.locSegment(1:end-1)*obj.rRotor).^2))';
            % segment radial length
            dR = obj.rRotor * (obj.locSegment(2:end)- ...
                 obj.locSegment(1:end-1))';
            
            % initial flow factors
            aSegment = 0.3 * ones(obj.nAnnulus, obj.nPsi);  
            apSegment = zeros(obj.nAnnulus, obj.nPsi);     
            
            for j = 1:obj.nIter
                % compute velocities for elements in annulus
                % axial velocity
                uRotor = obj.uInf*cos(obj.yawAngle)*(1-aSegment);
                % tang. velocity
                uTan = (1+apSegment).*(obj.Omega.*obj.rR*obj.rRotor ...
                    - obj.uInf*sin(obj.yawAngle)*cos(obj.psiSegment));
                % relative velocity
                uPer = sqrt(uRotor.^2 + uTan.^2);
                % inflow angles of each segment
                phiSegment = atan2(uRotor, uTan);
                % compute non-dim loads for elements in annulus
                [cAx, cAz, ClSegment, CdSegment, alphaSegment] ...
                    = obj.computeLoadsSegment(uRotor, uTan, twistAngle);
                % dimensionalize it but leave density out (cancels out)
                % axial force
                AxSegment = cAx.*0.5.*chord.*uPer.^2.*dR*obj.nBlades;
                % azimuthal force
                AzSegment = cAz.*0.5.*chord.*uPer.^2.*dR*obj.nBlades;
                % thrust coefficient
                CTSegment = AxSegment ./ (0.5*areaSegment*obj.uInf^2);     

                % log axial force each iteration
                obj.thrustIter(:,:,j) = AxSegment;
                
                % local solidity
                sigmaR = obj.nBlades*chord ./ (2*pi*obj.rR*obj.rRotor);
                fTotSegment = obj.calcPrandtlCorr(phiSegment);
                
                % compute a
                kappa   = sigmaR .*cAx ./ (4*fTotSegment.*sin(phiSegment).^2);
                kappap  = sigmaR .*cAz ./ (4*fTotSegment.*sin(phiSegment)...
                    .*cos(phiSegment));

                % if kappa > 2/3
                gamma1 = 2*fTotSegment.*kappa - (10/9-fTotSegment);
                gamma2 = 2*fTotSegment.*kappa - fTotSegment.* ...
                    (4/3-fTotSegment);
                gamma3 = 2*fTotSegment.*kappa - (25/9-2*fTotSegment);
                aip1Segment = (gamma1-sqrt(gamma2))./gamma3;
                % if kappa < 2/3
                aip1Segment(kappa <= 2/3) = kappa(kappa <=2/3) ./ ...
                    (1+kappa(kappa<=2/3));
                % if kappa > 1
                aip1Segment(kappa>1 & phiSegment<0) = kappa(kappa>1 & ...
                    phiSegment<0) ./ (kappa(kappa>1 & phiSegment<0)-1);
                
                % correct for skewed wake
                aip1Segment = obj.skewWakeCorr(aip1Segment);
                % limit updates
                aip1Segment(aip1Segment > 0.95) = 0.95;
                % compute a'
                apip1Segment = kappap ./ (1-kappap);
                
                aSegment = 0.75*aSegment + 0.25*aip1Segment;
                apSegment = 0.75*apSegment + 0.25*apip1Segment;
                
                % compute CT from a
                
                
                % finish iterating if error is below tolerance
                if all(abs(aSegment - aip1Segment) < obj.atol, "all") &&... 
                   all(abs(apSegment - apip1Segment) < obj.atol, "all")
                    break;
                end
            end
            
            % copy last value in col array to remove zeros
            obj.thrustIter(obj.thrustIter==0)= ...
                obj.thrustIter(find(obj.thrustIter,1,'last'));
            % compute circulation from uRotor and uTan
            obj.gamma = 0.5*sqrt(uRotor.^2+uTan.^2).*ClSegment.*chord/...
                (pi*obj.uInf^2/obj.Omega/obj.nBlades);
            
            % compute other parameters
            % torque coefficient
            CQSegment = AzSegment ./ (0.5*areaSegment*obj.uInf^2);
            % normal coefficient
            CNSegment = CQSegment ./ obj.rR;
            
            % store final result in corresponding property
            obj.areaAnnulus = areaSegment;
            obj.alpha(:, :) = alphaSegment;
            obj.phi(:, :)   = phiSegment;
            obj.aprime(:, :)= apSegment;
            obj.a(:, :)     = aSegment;
            obj.Ax(:, :)    = AxSegment;
            obj.Az(:, :)    = AzSegment;
            obj.Cl          = ClSegment;
            obj.Cd          = CdSegment;
            obj.CT(:, :)    = CTSegment;
            obj.CN(:, :)    = CNSegment;
            obj.Cq(:, :)    = CQSegment;
            obj.fTot(:, :)  = fTotSegment;
        end
    end
    
    methods(Static)
        function chordLength = computeChordLength(rR)
            % CHECKED
            % Determine the chord distribution for each blade segment
            chordLength = 3*(1-rR)+1;
        end
        
        function twistAngle = computeTwists(rR)
            % CHECKED
            % Determine the twist angle distribution for each blade segment
            twistAngle = deg2rad(14*(1-rR));
        end
    end
end

