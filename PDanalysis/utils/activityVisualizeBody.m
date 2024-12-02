function activityVisualizeBody(ID,fontsize)
%
% Description:
% This function generates a visualization of activity patterns for a given 
% participant based on their trunk sensor data. It provides a day-by-day 
% visualization where activities (e.g., walking, sitting/standing, lying, other) 
% are color-coded, and the corresponding acceleration data is overlaid. 
% The x-axis represents time within a 24-hour window, while each subplot 
% represents a single day of recorded data.
%
% Inputs:
% - ID: (string) Participant ID, e.g., '01'.
%       The function assumes that the corresponding processed data file 
%       (e.g., 'ID01_processed.mat') exists in the 'dataProcessed_Body' directory.
%
% - fontsize: (numeric) Font size for axis labels, tick labels, and other plot text.
%
% Outputs:
% - This function does not return variables. It produces a figure showing 
%   activity patterns with one subplot per day.

% Assumptions:
% - The input `.mat` file contains the following variables:
%   - `label`: Activity labels for each data point.
%   - `timeDatetime`: Timestamps corresponding to each data point.
%   - `datause`: Sensor data, e.g., acceleration.
%   - `timeNorm`: Normalized time used to calculate the sampling frequency.
%
% Activity Labels:
% - Activity labels are represented by specific numeric values:
%   - 1: Walking
%   - 2: Sit/Stand
%   - 3: Lying
%   - 4: Other

% Visualization:
% - Each subplot corresponds to one day of data.
% - The x-axis shows time in hours (e.g., '00:00', '02:00').
% - Activities are displayed as colored blocks, while acceleration data is shown 
%   as a black line.

% Example Usage:
% activityVisualizeBody('01', 12);
% This visualizes the activity data for participant 'ID01' with a font size of 12.

matName=strcat('ID',ID,'_processed');
load(fullfile('dataProcessed_Body',matName),"label",'timeDatetime','datause','timeNorm')
fs = 1 / (timeNorm(2) - timeNorm(1)); % Calculate sampling frequency (Hz)

time=timeDatetime;
time=time+hours(1);
monthNames = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};

% Find unique days in the timeDatetime
days = unique(dateshift(time, 'start', 'day'));
datestring=strings(size(days));
RangeIndex=false(length(time),3);
for i=1:length(days)
    datestring(i)=strcat( string(day(days(i))),{' '},monthNames{month(days(i))});

end

days=[days;days(end)+day(1)];
for i=1:(length(days)-1)
    RangeIndex(:,i)=(time >= days(i)) & (time <= days(i+1));
end

Nodata24h=24*60*60*fs;

figure
plotlabeluse=nan(Nodata24h,1);
plotdatause=nan(Nodata24h,1);

% Define the colors for each activity type
% activity_colors = {'#FFC374', '#7F9F80', '#124076', '#7F9F80', '#F9E897','#1A2130'};
activity_colors = {'#FFC374', '#7F9F80', '#124076', '#F9E897'};
for i =1:(length(days)-1)
    % for i =2:2
    subplot((length(days)-1),1,i)

    if i==1
        datatimeuse=time(RangeIndex(:,i));
        plotlabeluse((end-length(datatimeuse)+1):end)=label(RangeIndex(:,i));
        plotindex=[length(plotlabeluse)-length(datatimeuse)+1,length(plotlabeluse)];
        plotdatause((end-length(datatimeuse)+1):end)=datause(RangeIndex(:,i))/6;
        % Rnagestartl=find((RangeIndex(:,i)>=plotindex(1))==1,1,'first');
        % Rnageendl=find((RangeIndex(:,i)>=plotindex(1))==1,1,'last');
    elseif i==length(days)
        datatimeuse=time(RangeIndex(:,i));
        plotlabeluse(1:(length(datatimeuse)-1))=label(RangeIndex(:,i));
        plotdatause(1:(length(datatimeuse)-1))=datause(RangeIndex(:,i))/6;
        plotindex=[1,(length(datatimeuse)-1)];
    else
        datatimeuse=time(RangeIndex(:,i));
        plotlabeluse=label(RangeIndex(:,i));
        plotdatause=datause(RangeIndex(:,i))/6;
        plotindex=[1,length(label(RangeIndex(:,i)))];
    end


    % Calculate the differences between consecutive elements
    plotlabeluse2=plotlabeluse;
    plotlabeluse2(isnan(plotlabeluse2))=100;
    differences = diff(plotlabeluse2);

    % Find the indices where changes occur (where the difference is not zero)
    changeIndices = find(differences ~= 0);


    if i==1
        changeIndices=[changeIndices;length(plotlabeluse)];
    elseif i==length(days)
        changeIndices=[0;changeIndices;length(plotlabeluse)];
    else
        changeIndices=[0;changeIndices;(length(datatimeuse)-1)];
    end

    % y_pos = y_ticks(end-i+1)
    y_pos=1;
    plotindexstart=plotindex(1);
    Nodatapoint=0;
    hold on
    for j=1:(length(changeIndices)-1)
        plotrange= (changeIndices(j)+1:changeIndices(j+1));
        Nodatapoint=length(plotrange);
        % rectangle('Position', [plotindexstart, y_pos-0.4,Nodatapoint, 0.8], ...
        %     'FaceColor', activity_colors{unique(plotlabeluse(plotrange))}, 'EdgeColor', 'none');
        rectangle('Position', [plotindexstart, y_pos-0.4,Nodatapoint, 0.8], ...
            'FaceColor', activity_colors{unique(plotlabeluse(plotrange))}, 'EdgeColor', 'none');
        plotindexstart=plotindexstart+Nodatapoint;
        % colorsuse(j)=activity_colors{unique(plotlabeluse(plotrange))};
        % colorsuse(j)=activity_colors{unique(plotlabeluse(plotrange))};
    end


    if i==1
        % add a legend
        for j = 1:length(activity_colors)
            plot(NaN,NaN,'color', activity_colors{j}, 'LineWidth', 5);
        end

    end

    % plot acceleration
    plot(plotdatause+0.5,'k','LineWidth',1.5)
    hold off

    if i==1
        legend_labels = {'Walking', 'Sit/Stand', 'Lying', 'Other','Acc'};
        legend( legend_labels, 'Location', 'northoutside', 'Orientation', 'horizontal');
    end
    % legend(legend_labels)
    if i==(length(days)-1)
        % Set the x-axis ticks and labels
        x_ticks = 0:(2*60*60*fs):Nodata24h;
        x_labels = {'00:00', '02:00', '04:00', '06:00', '08:00', '10:00', '12:00', '14:00', '16:00', '18:00', '20:00', '22:00', '24:00'};

    else
        x_ticks = 0:(2*60*60*fs):Nodata24h;
        x_labels={};
    end
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);
    % Set the x-axis limits (start and end times in hours)
    xlim([0, Nodata24h]);
    ylim([0.5 1.5])
    % Set the y-axis ticks and labels
    y_ticks = 1;
    set(gca, 'YTick', y_ticks, 'YTickLabel', datestring(i), 'YGrid', 'on');
    datestring(i)
    set(gca, 'FontSize',fontsize);
    % ylim([])

    grid on
    ax=gca;
    ax.YGrid = 'off'; % This turns off y-axis grid lines, emphasizing the x-axis grid
    ax.XMinorGrid = 'on'; % Enable minor grid for x-axis
    % ax.MinorGridLineStyle = '--'; % Optional: Set the style of the minor grid lines
end

end