classdef BEMsolverNREL
    %BEMsolverNREL NREL's implementation of BEM
    
    % user-defined
    properties
        rho(1,1) double {mustBeNonnegative} = 1.225;      % density kg/m2
        uInf (1,1) double {mustBeNonnegative} = 10.0;     % m/s
        TSR (1,1) double {mustBeNonnegative} = 6.0;       % tip speed ratio
        rRotor (1,1) double {mustBeNonnegative} = 50.0;   % radius of rotor
        rTipRatio (1,1) double {mustBeNonnegative} = 1.0  % r_tip/R_rotor
        rRootRatio (1,1) double {mustBeNonnegative} = 0.2;% r_root/R_rotor
        nSegments (1,1) double {mustBeNonnegative} = 10;   % blade segments
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
    properties(Access = private)
        %twistDistribution (1,:) double; % twist angle at each blade segment
        %chordDistribution (1,:) double; % chord length of each blade seg.
        locSegment (1, :) double;       % radial length of each blade seg.
        alphaSegment (1, :) double;     % local angle-of-attack
        phiSegment (1, :) double;       % local inflow angle
    end
    
    properties(Access = protected)
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
            obj = createPolarSplines(obj);
            obj = computeSegments(obj);
            obj.Omega = obj.TSR * obj.uInf / obj.rRotor;
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
                segmentBounds = linspace(rStart, rEnd, obj.nSegments+1);
%                 lengthSegments = diff(segmentBounds);
                
            elseif obj.spacing == "cosine"
                segmentBounds = cosspace(rStart, rEnd, obj.nSegments+1);
%                 lengthSegments = diff(segmentBounds);
                
            else
                disp("unrecognized spacing, set to default")
                segmentBounds = linspace(rStart, rEnd, obj.nSegments+1);
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
        
        function [cl, cd] = computeLoadsSegment(obj, uRotor, uTan, ...
                              twistAngle, bladePitch, fCL, fCD)
            theta = deg2rad(bladePitch);            % blade pitch
            phi = atan2(uRotor, uTan);              % flow angle
            alpha = phi - twistAngle - theta;       % angle of attack
            
            if alpha > obj.amax
                cl = fCL(obj.amax); cd = fCD(obj.amax);
            elseif alpha < obj.amin
                cl = fCL(obj.amin); cd = fCD(obj.amin);
            else
                cl = fCL(alpha); cd = fCD(alpha);
            end
        end
        
        function solTotal = solveRotor(obj)
            % solve problem for each annulus
            solTotalArr = zeros(obj.nSegments, 7);
            askewTotal = zeros(obj.nSegments, obj.nPsi);
            for i = 1:obj.nSegments
                [sol, askew] = solveStreamtube(obj, ...
                    obj.locSegment(i), obj.locSegment(i+1), ...
                    obj.rRootRatio, obj.rTipRatio, obj.rRotor, ...
                    obj.uInf, obj.Omega, obj.nBlades);
                solTotalArr(i, :) = sol(:);
                askewTotal(i,:) = askew(:);
            
            solTotal = struct("rR", solTotalArr(:,1), "a", solTotalArr(:,2), ...
                    "aprime", solTotalArr(:,3), "askew", askewTotal, ...
                    "nAx", solTotalArr(:,4), "nAz", solTotalArr(:,5), ...
                    "fTot", solTotalArr(:,7), "CT", solTotalArr(:,6));
            end
        end
        
        function [sol, askew] = solveStreamtube( ...
                obj, rR1, rR2, rRoot, rTip, rRotor, uInf, Omega, nBlades)
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
            
            rR = (rR1 + rR2)/2;                 % center of blade segment
            twistAngle = obj.computeTwists(rR); % twist of blade seg.
            chord = obj.computeChordLength(rR);
            a = 0.3; aprime = 0;  % flow factors
            psiArr = linspace(0, 2*pi, obj.nPsi);
            for i = 1:obj.nIter
                
                % calculate velocity and loads 
                uRotor = uInf*(1-a);               % axial velocity
                uTan = (1+aprime)*Omega*rR*rRotor; % tangential velocity
%                 qDyn = 0.5*obj.rho*uPer^2;
                [cl, cd] = obj.computeLoadsSegment(uRotor, uTan, ...
                                   twistAngle, obj.bladePitch, ...
                                   obj.fCL, obj.fCD);
                
                % calculate inflow angle from velocities
                phi = atan2(uRotor, uTan);
                sigmaP = nBlades*chord/(2*pi*rR*rRotor);
                CT = 1+sigmaP*(1-a)^2*(cl*cos(phi)+cd*sin(phi))/ ... 
                     sin(phi)^2;
                
                % Tip and root correction
                fTot = obj.calcPrandtlTipCorr(rR, rRoot, rRotor, nBlades, ...
                       phi);
                 
                % compute a and aprime
                if CT > 0.96*fTot
                    a_ip1 = (18*fTot-20-3*sqrt(CT*(50-36*fTot)+12*fTot...
                            *(3*fTot-4)))/(36*fTot-50);
                else
                    a_ip1 = (-1+4*fTot*sin(phi)^2/ ...
                             (sigmaP*(cl*cos(phi)+cd*sin(phi))))^(-1);
                end
                
                aprime_ip1 = (-1+4*fTot*sin(phi)*cos(phi)/ ...
                             (sigmaP*(cl*sin(phi)-cd*cos(phi))))^(-1);
                
                if a_ip1 > 0.95
                    a_ip1 = 0.95;
                end
                
                if abs(a - a_ip1) < obj.atol && ... 
                   abs(aprime - aprime_ip1) < obj.atol
                    break;
                end
                
                % update flow factors
                a = 0.75 * a + 0.25 * a_ip1;
                aprime = 0.75 * aprime + 0.25 * aprime_ip1;
                
            end
            % apply skewing factor / correction for yaw
            askew = obj.skewWakeCorr(a, rR, obj.yawAngle, psiArr);
            sol = [rR, a, aprime, cl, cd, CT, fTot];
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
        
        function fTot = calcPrandtlTipCorr(rR, rRoot, rRotor, nBlades, phi)
            % CHECKED
            % calculate Prandtl Tip Corrections FACTORS!
            r = rR * rRotor;
            rHub = rRoot * rRotor;
            temp1 = -nBlades/2*(rRotor - r)/(r*sin(phi));
            fTip  = 2/pi*acos(exp(temp1));
            temp1 = -nBlades/2*(r - rHub)/(r*sin(phi));
            fRoot = 2/pi*acos(exp(temp1));
            fTot = fRoot*fTip;
            if fTot < 1e-4 
                fTot = 1e-4;  % avoide divide by zero or blow-up
            elseif isnan(fTot) == true
                fTot = 0;
            end
        end

        function a = calcGlauertCorr(CT)
            % CHECKED
            % computes the Glauert correction for heavily loaded rotors
            if CT < (2*sqrt(1.816)-1.816)
                a = 0.5 - sqrt(1-CT)/2;
            else
                a = 1+(CT-1.816)/(4*sqrt(1.816)-4);
            end
        end
        
        function aSkew = skewWakeCorr(a, rR, gamma, psi)
            chi = gamma*(1+0.6*a);
            aSkew = a*(1+15*pi/32*rR*tan(chi/2).*cos(psi));
        end
    end
end

