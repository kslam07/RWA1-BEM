function [out, R1, R2, R3, R4]=enthalpy(solver)
    p=101325;
    rho=1.225;
    pt=p+0.5*rho*solver.uInf^2;
    v2=solver.uInf*(1-mean(solver.a,2));
    p2=pt-0.5*rho*v2.^2;

    pressureJump=sum(solver.Ax*solver.rho,2)./solver.areaAnnulus;
    v4=solver.uInf*(1-2*mean(solver.a,2));
    p3=p2-pressureJump;
    pt3=p3+0.5*rho*v2.^2;
    p4=pt3-0.5*rho*v4.^2;

    % compute radius at stations
    R1 = sqrt(sum(solver.areaAnnulus.*v2/solver.uInf)/pi);
    R2 = solver.rRotor;
    R3 = solver.rRotor;
    R4 = sqrt(sum(solver.areaAnnulus.*v2./v4)/pi);
    
    enthalpy1=ones(1,solver.nAnnulus)*(p/rho+solver.uInf^2/2);
    enthalpy2=transpose(mean(ones(1,solver.nAnnulus).*(p2/rho+v2.^2/2),2));
    enthalpy3=transpose(mean(ones(1,solver.nAnnulus).*(p3/rho+v2.^2/2),2));
    enthalpy4=transpose(mean(ones(1,solver.nAnnulus).*(p4/rho+v4.^2/2),2));
    enthalpy3mean=ones(1,solver.nAnnulus).*mean(mean(enthalpy3,2));
    out=cat(1,enthalpy1, enthalpy2, enthalpy3, enthalpy4, enthalpy3mean);
end 