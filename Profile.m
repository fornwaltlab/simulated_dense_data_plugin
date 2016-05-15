classdef Profile < hgsetget
    % Profile - GUI Element for manipulating a plot
    %
    %   This class creates a plot object that can be changed by dragging
    %   control points around. This isn't a contour, but rather a simple
    %   profile. It was designed explicitly for DNSCreator generation of
    %   radial and twist profiles.

    properties
        XData                       % X Coordinates of the reference points
        YData                       % Y Coordinates of the reference points

        Parent                      % Axes which will contain the plot
        MarkerSize  = 20;           % Size of the Markers to use
        LineWidth   = 2;            % Width of the line connecting points
        Color       = 'r';          % Color of the profile curve

        XLimits     = [-inf inf];   % Constrains on X position
        YLimits     = [-inf inf];   % Constrains on Y position
    end

    properties (Hidden)
        hPoints                     % Handle to the control points
        hLine                       % Handle to the line object
        hGroup                      % Handle to the whole ensemble
    end

    properties (Dependent)
        Points                      % Returns the coordinates of the profile
    end

    events
        EVT_Changed                 % Event to be triggered on changes
    end

    methods

        %--- Getter Functions ---%
        function res = get.Points(self)
            res = [self.XData(:), self.YData(:)];
        end

        %--- Getter Functions ---%
        function set.Points(self, val)
            self.XData = val(:,1);
            self.YData = val(:,2);

            refresh(self);
        end

        %--- Constructor ---%
        function self = Profile(xdata, ydata, parent, varargin)
            % Profile - Constructor for profile object
            %
            % USAGE:
            %   self = Profile(xdata, ydata, parent)
            %
            % INPUTS:
            %   xdata:  [M x 1], X Coordinates of control points
            %   ydata:  [M x 1], Y Coordinates of control points
            %   parent: Handle, Axes in which to place all elements
            %
            % OUTPUTS:
            %   self:   Object, Handle to the graphics object to manipulate
            %           display and return values

            self.XData = xdata;
            self.YData = ydata;

            if ~exist('parent', 'var'); parent = gca; end

            % If any addition inputs were provided, go ahead and set them
            if numel(varargin); set(self,varargin{:}); end

            self.Parent = parent;

            % Initialize the actual graphics object
            self.init()

            % Ensure that the plots are up-to-date
            self.refresh()
        end

        function clickLine(self)
            % clickLine - Callback for clicking on the line to add point
            %
            % USAGE:
            %   self.clickLine()

            cp = get(gca, 'currentpoint');

            % Add the point ta the end
            self.XData(end+1) = cp(1,1);
            self.YData(end+1) = cp(1,2);

            % Sort to ensure monotonicity
            [self.XData, sortind] = sort(self.XData);
            self.YData = sort(self.YData);
            [~,ind] = max(sortind);

            % Now allow the user to start dragging without clicking again
            self.dragPoint(ind);

            % Update display (of course)
            refresh(self);
        end

        function clickPoint(self)
            % clickPoint - Callback for clicking a control point

            % Get the current point in the axes
            cp = get(gca, 'currentpoint');

            % Determine the closest neighbor
            [~,ind] = min(abs(self.XData - cp(1)));

            switch lower(get(gcf, 'selectiontype'))
                case 'normal'                       % Drag the point
                    self.dragPoint(ind);
                case 'alt'                          % Remove the point
                    self.deletePoint(ind);
                    notify(self, 'EVT_Changed');
                case 'open'                         % Specify Location
                    self.specifyPoint(ind);
                    notify(self, 'EVT_Changed');
            end
        end

        function deletePoint(self, ind)
            % deletePoint - Removes the specified point from the profile
            %
            % USAGE:
            %   self.deletePoint(ind)
            %
            % INPUTS:
            %   ind:    Integer, Index of the point to remove

            % Don't allow the user to remove end points as there must be at
            % least two points
            if ind == numel(self.XData) || ind == 1
                return;
            end
            self.XData(ind) = [];
            self.YData(ind) = [];

            % Update display
            refresh(self);
        end

        function drag(self, validx, validy, ind)
            % drag - Point drag callback executed when mouse moves
            %
            % USAGE:
            %   self.drag(validx, validy, ind)
            %
            % INPUTS:
            %   validx: Functionhandle, function used to validate the new
            %           X position
            %   validy: Functionhandle, function used to validate the new
            %           Y position
            %   ind:    Integer, index of the point that we are dragging

            cp = get(gca, 'currentpoint');
            self.XData(ind) = validx(cp(1,1));
            self.YData(ind) = validy(cp(1,2));

            % Update the display to reflect the drag
            refresh(self)
        end

        function dragPoint(self, ind)
            % dragPoint - Allows the user to drag the specified point
            %
            % USAGE:
            %   self.dragPoint(ind)
            %
            % INPUTS:
            %   ind:    Integer, Index of the point to drag around

            % Create some constraints on where the point can be dragged to
            if ind == 1
                validx = @(x)min(self.XData);
            elseif ind == numel(self.XData)
                validx = @(x)max(self.XData);
            else
                % Points cannot be reordered
                validx = @(x)self.constrain(x, self.XData(ind-1),...
                                               self.XData(ind+1));
            end

            validy = @(y)min(max(self.YLimits(1), y), self.YLimits(2));

            fig = ancestor(self.Parent, 'figure');

            % Remember the old figure properties so we can reset them
            opts.WindowButtonMotionFcn = get(fig, 'WindowButtonMotionFcn');
            opts.WindowButtonUpFcn = get(fig, 'WindowButtonUpFcn');

            % Set Motion Callbacks for drag event
            set(fig,'WindowButtonUpFcn', @(varargin)self.undrag(fig, opts));
            set(fig,'WindowButtonMotionFcn', @(s,e)self.drag(validx,validy,ind));
        end

        function refresh(self)
            % refresh - Redraws the profile according to current points
            %
            % USAGE:
            %   self.refresh()
            %

            set(self.hPoints, 'xdata', self.XData, 'ydata', self.YData)

            % Interpolate the line
            xx = linspace(self.XData(1),self.XData(end),1000);
            pp = spline(self.XData, self.YData);
            yy = ppval(pp, xx);

            set(self.hLine, 'xdata', xx, 'ydata', yy)
        end

        function specifyPoint(self, ind)
            % specifyPoint - Allows the user to define EXACT point location
            %
            % USAGE:
            %   self.specifyPoint(ind)
            %
            % INPUTS:
            %   ind:    Integer, Index of the point that we are going to
            %           move to the specified location

            % If the user selects an end point only allow them to move it
            % up and down
            if ind == 1 || ind == numel(self.XData)
                prompt = {'Y:'};
                def = {num2str(self.YData(ind))};
            else
                prompt = {'X:','Y:'};
                def = {num2str(self.XData(ind)), num2str(self.YData(ind))};
            end
            title = 'Specify Coordinate';
            lines = 1;

            % Prompt the user for the new point location
            answer = inputdlg(prompt, title, lines, def);

            % If the user hits cancel, then just skip the rest
            if isempty(answer); return; end

            % Ensure the points are valid
            if numel(answer) == 2
                x = eval(answer{1});

                if x < self.XLimits(1) || x > self.XLimits(2)
                    errordlg(sprintf('X Values must be between %d and %d', ...
                        self.XLimits(1), self.XLimits(2)));
                    return;
                end

                self.XData(ind) = x;
            end

            y = eval(answer{end});
            if y < self.YLimits(1) || y > self.YLimits(2)
                errordlg(sprintf('Y Values must be between %d and %d', ...
                    self.YLimits(1), self.YLimits(2)));
                return;
            end

            self.YData(ind) = y;

            % Make sure that we have a monotonic function by sorting the
            % results
            [self.XData, sortind] = sort(self.XData);
            self.YData = self.YData(sortind);

            % Update the display to reflect changes
            refresh(self);
        end

        function undrag(self, hfig, opts)
            % undrag - Helper function for stopping drag events
            %
            % USAGE:
            %   self.undrag(fig, opts)
            %
            % INPUTS:
            %   fig:    Handle, handle to the figure controlling the
            %           current drag event
            %   opts:   Properties and values to set for this figure

            set(hfig, opts)

            % Trigger the change event so listeners can update
            notify(self, 'EVT_Changed');
        end
    end

    methods (Static)
        function x = constrain(x, xlow, xhigh)
            % constrain - Helper function for constraining coordinates
            if x <= xlow
                x = xlow + eps;
            elseif x >= xhigh
                x = xhigh - eps;
            end
        end
    end

    methods (Access = 'private')
        function init(self)
            % init - Initialize all GUI components
            %
            %   This function should only be called on object creation.
            %
            % USAGE:
            %   self.init()

            self.hGroup = hggroup('Parent', self.Parent);

            self.hLine = line(NaN, NaN, 'Parent',       self.hGroup,...
                                        'linestyle',    '-',...
                                        'marker',       'none',...
                                        'color',        self.Color,...
                                        'linewidth',    self.LineWidth);

            self.hPoints = line(NaN, NaN, 'Parent',     self.hGroup,...
                                        'linestyle',    'none',...
                                        'marker',       '.',...
                                        'markersize',   20,...
                                        'color',        self.Color,...
                                        'Markersize',   self.MarkerSize);

            set(self.hPoints, 'buttondownfcn', @(s,e)self.clickPoint());
            set(self.hLine, 'buttondownfcn', @(s,e)self.clickLine());
        end
    end
end
