%% Original version of the bump measurement
% Detect if there was a bump at the end of the simulation
% Measure activity bump width as width at 50% of maximum
function bump=measureBump1(lastActivity, fixedThreshold)
    if fixedThreshold
        lastActivity=lastActivity>1e-6;
    else
        activityThreshold = (max(lastActivity) - min(lastActivity)) / 2;
        lastActivity=lastActivity>activityThreshold;
    end
    bumpMargins = diff([0 lastActivity 0]);
    firstBumpMargins = find(bumpMargins~=0, 2);
    if size(firstBumpMargins, 2) < 2 % No bump was found
        bump.bumpWidth = 0;
        bump.bumpLoc   = 0;
    else
        bumpStart = firstBumpMargins(1);
        bumpEnd   = firstBumpMargins(2);
        % Calculate the bump width in rads
        bump.bumpWidth = (((bumpEnd-1) - bumpStart) / 8) * 2 * pi;
        bump.bumpLoc   = (((bumpEnd-1) + bumpStart) / 2);
    end

end