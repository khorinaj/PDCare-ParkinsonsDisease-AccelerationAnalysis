function [validStepTimes, validAccHS] = filterSteps(timeHS, accHS)
    % Calculate the step times (difference between consecutive heel strikes)
    stepTimes = diff(timeHS);

    % Identify valid steps (those less than or equal to 1.5 seconds)
    validSteps = stepTimes <= 1.5;

    % Initialize arrays to store valid step times and valid HS accelerations
    validStepTimes = [];
    validAccHS = [];

    % Find indices of rest periods
    restIndices = find(~validSteps);

    % Add boundary conditions to handle the start and end of the data
    restIndices = [0; restIndices; length(stepTimes)];

    % Loop through segments between rest periods
    for i = 1:length(restIndices)-1
        % Define the segment of interest
        segment = restIndices(i)+1 : restIndices(i+1);

        % Check if the segment has at least 5 valid steps
        if sum(validSteps(segment)) >= 5
            % Append valid steps in this segment to the result arrays
            validStepTimes = [validStepTimes; stepTimes(segment(validSteps(segment)))];
            validAccHS = [validAccHS; accHS(segment(validSteps(segment)))];
        end
    end
end
