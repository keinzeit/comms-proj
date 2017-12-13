function targetPulse = normalizePulse(targetPulse,sourcePulse)
% function targetPulse = normalizePulse(targetPulse,sourcePulse)
% Function normalizes the energy in the targetPulse to the energy in the
% sourcePulse. If the energy in targetPulse is given by
% sum(targetPulse.^2) and that in sourcePulse is given by
% sum(sourcePulse.^2), then the function multiplies targetPulse by a factor
% A such that sum((A*targetPulse.^2)) = sum(sourcePulse.^2).

tEnergy = sum(targetPulse.^2);
sEnergy = sum(sourcePulse.^2);

A = sqrt(sEnergy/tEnergy);

targetPulse = A*targetPulse;

return
