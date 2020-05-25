function break_axis(varargin)
%BREAK_AXIS Inserts a break on the selected axis.
%   BREAK_AXIS() applies a break to the y-axis of the current axis with 
%       default height and location.
%   BREAK_AXIS('handle', ax) applies a break to the axis ax.
%   BREAK_AXIS('length', L) applies a break with length L (L is a fraction
%       between ticks).
%   BREAK_AXIS('axis', 'x') applies a break to the x-axis.
%   BREAK_AXIS('bottom_reset', X) resets the bottom tick label to X (X 
%      must be a string).
%   BREAK_AXIS('position', P) breaks the axes at position P.

%% Parse the inputs
% Default values
ax = gca;
extent = 0.5;
width = get(ax, 'ticklength');
width = width(1)*2;
bottom_reset = '';
xy = 'y';
position = Inf;

% Overwrite with what was fed in
if nargin > 0
    for i=1:2:length(varargin)
        if strcmp(varargin(i), 'handle')
            %ax = cell2mat(varargin(i+1));
            ax = varargin(i+1);
            ax = ax{1};
        end
        if strcmp(varargin(i), 'length')
            extent = cell2mat(varargin(i+1));
        end
        if strcmp(varargin(i), 'axis')
            xy = cell2mat(varargin(i+1));
        end
        if strcmp(varargin(i), 'bottom_reset')
            bottom_reset = cell2mat(varargin(i+1));
        end
        if strcmp(varargin(i), 'position')
            position = cell2mat(varargin(i+1));
        end
    end
end
       
%% Do the things

% Get figure handle
f = get(ax, 'parent');

% Get number of ticks
n = length(get(ax, strcat(xy, 'ticklabel'))) - 1;
pos = get(ax, 'position');

if strcmp(xy, 'x')
    
    % Get the spacing between ticks
    spacing = pos(3)/n;
    %spacing = extent;
    
    % Find the center along the axis
    cent = pos(2) + 0.5*spacing;
    
    % Reset if a position was given
    if position < Inf
        lims = get(gca, 'xlim');
        cent = pos(1) + pos(3) * (position - lims(1))/(diff(lims));
    end
    
    % define the extents of the bar
    X = [cent - extent*spacing/2; 
         cent + extent*spacing/2];
    Y1 = [pos(2) - 0.005; 
          pos(2) + 0.005];
    Y2 = [pos(2) - width/2; 
          pos(2) + width/2];
    
    % Make the bar and ticks
    h1 = annotation(f, 'rectangle', [X(1), Y1(1), diff(X), diff(Y1)]);
    annotation(f, 'line', [X(1) X(1)], [Y2(1) Y2(2)]);
    annotation(f, 'line', [X(2) X(2)], [Y2(1) Y2(2)]);
    
elseif strcmp(xy, 'y')
    % Get the spacing between ticks
    spacing = pos(4)/n;
    
    % Find the center along the axis
    cent = pos(2) + 0.5*spacing;
    
    % Reset if a position was given
    if position < Inf
        lims = get(gca, 'ylim');
        cent = pos(2) + pos(4) * (position - lims(1))/(diff(lims));
    end
    
    % define the extents of the bar
    X1 = [pos(1) - 0.005; 
          pos(1) + 0.005];
    X2 = [pos(1) - width/2; 
          pos(1) + width/2];
    Y = [cent - extent*spacing/2; 
         cent + extent*spacing/2];
     
    % Make the bar and ticks
    h1 = annotation(f, 'rectangle', [X1(1), Y(1), diff(X1), diff(Y)]);
    annotation(f, 'line', [X2(1) X2(2)], [Y(1) Y(1)]);
    annotation(f, 'line', [X2(1) X2(2)], [Y(2) Y(2)]);
    
end

% Reset the color
h1.FaceColor = 'w';
h1.Color = 'w';

% Reset the bottom value
%if strcmp(bottom_reset, '')
%    list = get(ax, [xy 'ticklabel']);
%    list(1) = {bottom_reset};
%    set(gca, [xy 'ticklabel'], list);
%end

