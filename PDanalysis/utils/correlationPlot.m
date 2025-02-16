function correlationPlot(values,pvalues,savepath,varargin)

for i = 1:numel(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i}, 'labels')
        % Look for 'method' in varargin and set the method variable
        if i < numel(varargin)
            xlabel = varargin{i + 1};
            ylabel = varargin{i + 2};
            break
        end

    end
end
% Create the heatmap
fig = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full-screen figure
h = heatmap(strrep(cellstr(xlabel), '_', '\_'), strrep(cellstr(ylabel), '_', '\_'), ...
    values, 'MissingDataColor', 'w', 'GridVisible', 'on', 'CellLabelColor', 'None');

% Set the colormap
colormap(h, cool);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18);
h.ColorLimits = [-1, 1];
% Rotate X-axis labels manually
ax = struct(h); % Access the underlying axes
ax.XAxis.TickLabelRotation = 45; % Rotate the X-axis labels by 45 degrees

% Retrieve the heatmap's underlying axes
ax = gca;

% Set the figure's SizeChangedFcn to dynamically update annotations
fig.SizeChangedFcn = @updateAnnotations;

% Initial call to ensure annotations are added at the start
updateAnnotations();

% Callback function to dynamically adjust annotations
    function updateAnnotations(~, ~)
        % Clear existing annotations
        delete(findall(gcf, 'Tag', 'HighlightRectangle'));

        % Get the heatmap's position
        heatmapPosition = ax.Position; % [x, y, width, height]

        % Get the number of rows and columns
        numRows = size(values, 1);
        numCols = size(values, 2);

        % Identify cells where pvalues < 0.05
        [row, col] = find(pvalues < 0.05);

        % Loop to add rectangles
        for i = 1:length(row)
            % Calculate normalized coordinates for each rectangle
            xStart = heatmapPosition(1) + heatmapPosition(3) * (col(i) - 1) / numCols;
            yStart = heatmapPosition(2) + heatmapPosition(4) * (numRows - row(i)) / numRows; % Invert Y-axis
            width = heatmapPosition(3) / numCols;
            height = heatmapPosition(4) / numRows;

            % Draw a rectangle annotation
            annotation('rectangle', [xStart, yStart, width, height], ...
                'Color', 'k', 'LineWidth', 2, 'Tag', 'HighlightRectangle');
        end
    end


% Save the figure as a PNG file
saveas(gcf, savepath);
fprintf('Saved correlation Map in %s.\n', savepath);

end
