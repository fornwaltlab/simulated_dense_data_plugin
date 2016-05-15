classdef DNSCreator < hgsetget
    % DNSCreator - Graphically generates DNS files
    %
    %   This program supplies the user with an interface to changes various
    %   factors which affect DENSE image quality and analysis. Using this
    %   software, it is possible to generate a .dns file which can be passed
    %   directly to DENSEanalysis for strain calculation.

    properties (Dependent, SetObservable)
        DENC            % Displacement encoding values (cyc / mm)
        Frames          % Number of cardiac phases
        InnerRadius     % Endocardial radius (mm)
        OuterRadius     % Epicardial radius (mm)
        PixelHeight     % Vertical Pixel Spacing (mm / px)
        PixelSpacing    % DICOM-style Pixel Spacing field
        PixelWidth      % Horizontal Pixel Spacing (mm / px)
        RadialProfile   % Radial Contraction profile
        TwistProfile    % Twist profile from endo to epi
    end

    properties (Hidden)
        DNS                     % Variable for storing DNS structure
        Frame           = 1;    % Frame for preview display
        Handles                 % Handles to all graphics objects
        Padding         = 0.15; % Padding around epicardium in images
        Timer                   % Timer object for displaying movies
        Listen                  % Listener for destruction
        Data                    % DENSEdata object
        Loading                 % Boolean specifying whether to ignore
                                % preview or not
    end

    events
        StateChanged    % Event fired when anything changes
    end

    methods
        %--- Property Getters ---%
        function h = get.RadialProfile(obj); h = obj.Handles.radprof.Points; end
        function h = get.TwistProfile(obj);  h = obj.Handles.tprof.Points; end
        function h = get.InnerRadius(obj);   h = obj.val('InnerRadius'); end
        function h = get.OuterRadius(obj);   h = obj.val('OuterRadius'); end
        function h = get.Frames(obj);        h = obj.val('Frames'); end
        function h = get.PixelWidth(obj);    h = obj.val('PixelWidth'); end
        function h = get.PixelHeight(obj);   h = obj.val('PixelHeight'); end
        function h = get.DENC(obj);          h = obj.val('DENC'); end

        function h = get.PixelSpacing(obj)
            h = [obj.PixelHeight, obj.PixelWidth];
        end

        %--- Property Setters ---%
        function set.RadialProfile(self, val)
            self.Handles.radprof.Points = val;
            notify(self, 'StateChanged');
        end

        function set.TwistProfile(self, val)
            self.Handles.tprof.Points = val;
            notify(self, 'StateChanged');
        end

        function set.InnerRadius(self, val);   self.setVal('InnerRadius', val); end
        function set.OuterRadius(self, val);   self.setVal('OuterRadius', val); end
        function set.Frames(self, val);        self.setVal('Frames', val); end
        function set.PixelWidth(self, val);    self.setVal('PixelWidth', val); end
        function set.PixelHeight(self, val);   self.setVal('PixelHeight', val); end
        function set.DENC(self, val);          self.setVal('DENC', val); end

        function set.PixelSpacing(self, val)
            set(self, 'PixelHeight', val(1), 'PixelWidth', val(2))
        end

        %--- Object Constructor ---%
        function self = DNSCreator(varargin)
            % DNSCreator - Constructor for DNS-generating GUI
            %
            % USAGE:
            %   self = DNSCreator()
            %
            % OUTPUTS:
            %   self:   Object, Handle to the GUI which can be used to
            %           manipulate the data or graphical interface

            % Initialize all GUI components
            init(self);

            % Create timer object for movies
            self.Timer = timer('ExecutionMode', 'FixedRate',...
                               'ErrorFcn',      '',...
                               'Period',        0.1,...
                               'TimerFcn',      @(varargin)self.increment);

            % Add a listener to ensure that the timer gets deleted
            self.Listen = addlistener(self, 'ObjectBeingDestroyed',...
                                            @(varargin)delete(self.Timer));

            self.Loading = true;
            set(self, varargin{:})
            self.Loading = false;

            % Go ahead and preview the basic configuration
            preview(self);
        end

        function decrement(self)
            % decrement - Display the previous frame
            %
            % USAGE:
            %   self.decrement()

            self.Frame = mod(self.Frame - 2, self.Frames) + 1;
            self.gotoFrame(self.Frame);
        end

        function increment(self)
            % increment - Display the previous frame
            %
            % USAGE:
            %   self.increment()

            self.Frame = mod(self.Frame, self.Frames) + 1;
            self.gotoFrame(self.Frame);
        end

        function keypress(self, evnt)
            % keypress - Callback for processing all keypress events
            %
            % USAGE:
            %   self.keypress(evnt)
            %
            % INPUTS:
            %   evnt:   EventData, Object containing information about the
            %           key press event including the key that was pressed.

            switch lower(evnt.Key)
                case {'n', 'rightarrow'}  % Displays the next frame
                    self.increment();
                case {'b', 'leftarrow'}   % Displays the previous frame
                    self.decrement();
                case 'space'              % Stops and starts the movie
                    if strcmpi(self.Timer.Running, 'on')
                        stop(self.Timer);
                    else
                        start(self.Timer);
                    end
            end
        end

        function preview(self)
            % preview - Processes and displays the current DENSE parameters
            %
            %   This function is used to generate DENSE data using the current
            %   parameters and then display the phase images in the GUI.
            %
            % USAGE:
            %   self.preview()

            if self.Loading
                return;
            end

            self.DNS = generatedns(struct(self));

            % Ensure that the current frame isn't out of range
            self.Frame = min(self.Frame, self.DNS.dns.Number);

            % Now display it
            self.gotoFrame(self.Frame);
        end

        function save(self, filename)
            % save - Saves the current configuration to a DNS file
            %
            %   This function generates the DENSE images and contours based
            %   upon the current combination of parameters. The result is then
            %   saved to a DNS file which can then be loaded directly into
            %   DENSEanalysis.
            %
            % USAGE:
            %   self.save(filename)
            %
            % INPUTS:
            %   filename:   String, (Optional) Path to the .dns file to save
            %               the results in. If not supplied, a dialog box will
            %               be provided to allow the user to select a file.

            % If Data is defined, just load it into that
            if isvalid(self.Data) && ~isequal(self.Data, [])
                self.Data.load(@(s,p)deal(self.DNS, '', ''));
                return;
            end

            if ~exist('filename', 'var')
                [fname, pname] = uiputfile({'*.dns','DENSEanalysis File'});

                % Determine if the user hit cancel
                if isequal(fname, 0) || isequal(pname, 0); return; end

                filename = fullfile(pname, fname);
            end

            dns = self.DNS; %#ok
            save(filename, '-struct', 'dns')
        end

        function res = struct(self)
            % struct - Converts the object to a structure representation
            %
            % USAGE:
            %   S = struct(self)
            %
            % OUTPUTS:
            %   S:  Structure, Fields corresponding to fields within the
            %       object

            res = struct('InnerRadius',     self.InnerRadius,...
                         'OuterRadius',     self.OuterRadius,...
                         'RadialProfile',   self.RadialProfile,...
                         'TwistProfile',    self.TwistProfile,...
                         'PixelSpacing',    self.PixelSpacing,...
                         'DENC',            self.DENC,...
                         'Frames',          self.Frames,...
                         'Padding',         self.Padding);
        end
    end

    methods (Access = 'private')
        function gotoFrame(self, ind)
            % gotoFrame - Helper function for displaying data
            %
            %   This function performs the refresh of the imaging data based
            %   upon the currently specified frame
            %
            % USAGE:
            %   self.gotoFrame(ind)
            %
            % INPUTS:
            %   ind:    Integer, Frame to display in the current figure

            set(self.Handles.previewxim, 'CData', self.DNS.img{2}(:,:,ind));
            set(self.Handles.previewyim, 'CData', self.DNS.img{3}(:,:,ind));

            % Make sure that the axes are configured properly
            axis([self.Handles.previewx, self.Handles.previewy], 'tight');

            % Update the Frame Number label in the lower lefthand corner
            set(self.Handles.ftext, 'string', {'',['    ',num2str(ind)]});
        end

        function init(self)
            % init - Initialize all GUI components
            %
            %   This function contains all GUI-generating code and should only
            %   be called on object construction.
            %
            % USAGE:
            %   self.init()

            h.fig = figure('Toolbar',       'none',...
                           'Menubar',       'none',...
                           'Name',          'DNSCreator',...
                           'NumberTitle',   'off',...
                           'DeleteFcn',     @(varargin)delete(self));

            colormap gray;

            width   = 1200;
            height  = 400;

            mpos = get(0, 'monitorposition'); mpos = mpos(1,:);
            pos = [mpos(1) + ((mpos(3) - mpos(1))/2) - (width / 2),...
                   mpos(2) + ((mpos(4) - mpos(2))/2) - (height / 2),...
                   width, height];

            set(h.fig, 'position', pos);

            % Create some panels to keep shit organized better
            h.settings  = uipanel('Parent',     h.fig,...
                                  'Position',   [0.01 0.11 0.33 0.83]);
            h.cardiac   = uipanel('Parent',     h.settings,...
                                  'Title',      'Cardiac Parameters',...
                                  'Position',   [0.03 0.6 0.44 0.37]);
            h.imaging   = uipanel('Parent', h.settings,...
                                  'Position',   [0.53 0.6 0.44 0.37],...
                                  'Title',  'Imaging Parameters');

            % Axes for displaying the X-encoded phase data
            h.previewx = axes('parent',     h.fig,...
                               'box',       'on',...
                               'position',  [0.35 0.05 0.31 0.9]);

            h.previewxim = imagesc([], 'Parent', h.previewx);
            set(h.previewx, 'Xtick', [], 'ytick', [], 'ydir', 'reverse')

            % Set the CLims to the full dynamic range of a 12-bit dataset
            set(h.previewx, 'clim', [0, 4096])

            hold(h.previewx, 'on')

            % Frame number display
            h.ftext(1) = text(1, 1, '');

            % Axes for displaying the Y-encoded phase data
            h.previewy = axes('parent',     h.fig,...
                               'box',       'on',...
                               'position',  [0.68 0.05 0.31 0.9]);

            h.previewyim = imagesc([], 'Parent', h.previewy);
            set(h.previewy, 'Xtick', [], 'ytick', [], 'ydir', 'reverse')

            % Set the CLims to the full dynamic range of a 12-bit dataset
            set(h.previewy, 'clim', [0, 4096])

            hold(h.previewy, 'on')

            axis([h.previewx, h.previewy], 'square')

            % Frame number display
            h.ftext(2) = text(1, 1, '');

            % Ensure that the frame numbers don't look terrible
            set(h.ftext, 'FontWeight',          'bold',...
                         'VerticalAlignment',   'top');

            % Create a UIFLOWCONTAINER for buttons to minimize manual coding
            h.buttonfc = uiflowcontainer('v0',...
                                         'Parent',   h.fig,...
                                         'Position', [0.01 0.01 0.33 0.09]);

            % Use a margin of 5 pixels around buttons
            set(h.buttonfc, 'margin', 5)

            h.reset     = uicontrol('Parent',   h.buttonfc,...
                                    'String',   'Preview',...
                                    'Callback', @(varargin)self.preview());

            h.generate  = uicontrol('Parent',   h.buttonfc,...
                                    'Callback', @(varargin)self.save(), ...
                                    'String',   'Generate');

            % Create a UIGRIDCONTAINER for buttons to minimize manual coding
            h.cardiacfc = uigridcontainer('v0', 'Parent', h.cardiac);
            set(h.cardiacfc, 'GridSize', [3 2]);
            set(h.cardiacfc, 'Margin', 10);
            set(h.cardiacfc, 'HorizontalWeight', [2 1])

            % Create cardiac parameter items
            labels  = {'Endo Radius (mm)'
                       'Epi Radius (mm)'
                       'Phases'};

            tags    = {'InnerRadius'
                       'OuterRadius'
                       'Frames'};

            % Default values
            values  = [30 40 17];

            for k = 1:numel(labels)
                uicontrol('Parent', h.cardiacfc, ...
                          'Style',  'text',...
                          'String', labels{k});
                h.(tags{k}) = uicontrol('Parent', h.cardiacfc, ...
                                        'Style',  'edit',...
                                        'String', values(k));
            end

            h.imagingfc = uigridcontainer('v0', 'Parent', h.imaging);
            set(h.imagingfc, 'GridSize', [3 2]);
            set(h.imagingfc, 'Margin', 10);
            set(h.imagingfc, 'HorizontalWeight', [2 1])

            % Create cardiac parameter items
            labels  = {'Pixel Height (mm)'
                       'Pixel Width (mm)'
                       'Disp Enc. (cyc/mm)'};

            tags    = {'PixelHeight'
                       'PixelWidth'
                       'DENC'};

            % Default values
            values  = [1 1 0.1];

            for k = 1:numel(labels)
                uicontrol('Parent', h.imagingfc, ...
                          'Style',  'text',...
                          'String', labels{k});
                h.(tags{k}) = uicontrol('Parent', h.imagingfc, ...
                          'Style',  'edit',...
                          'String', values(k));
            end

            % Now create the axes with the radial and twist profiles
            h.rax = axes('Parent',      h.settings,...
                         'Position',    [0.03 0.1 0.44 0.37],...
                         'Xlimmode',    'manual');
            title(h.rax, 'Radial Profile (mm)');

            set(h.rax, 'XTick', [0.1 0.9], 'XTickLabel', {'ENDO', 'EPI'})

            % Creates an adjustable profile using PROFILE (external class)
            h.radprof = plugins.simulation.Profile([0 0.5 1], [10 3 0], h.rax);
            h.radprof.XLimits = [0 1];
            h.radprof.YLimits = [0 15];
            h.tlisten = addlistener(h.radprof, 'EVT_Changed', ...
                                    @(varargin)self.preview());

            h.tax = axes('Parent',      h.settings,...
                         'XLimMode',    'manual',...
                         'Position',    [0.53 0.1 0.44 0.37]);

            title(h.tax, 'Twist Profile (rad)');
            set(h.tax, 'XTick', [0.1 0.9], 'XTickLabel', {'ENDO', 'EPI'})
            set(h.tax, 'ylim', [-pi/4, pi/4])

            h.tprof = plugins.simulation.Profile([0 0.5 1], [pi/6 pi/12 0], h.tax);
            h.tprof.XLimits = [0 1];
            h.tprof.YLimits = [-pi/4 pi/4];

            h.tlisten = addlistener(h.tprof, 'EVT_Changed',...
                                    @(varargin)self.preview());

            % Store all graphics handles in the object
            self.Handles = h;

            % Ensure that changes to any edit box result in a refresh
            objs = findobj(h.fig, 'type', 'uicontrol', 'style', 'edit');
            set(objs, 'callback', @(varargin)self.preview())

            % Make sure that we set the keypress for all objects
            f = findobj(h.fig, '-property', 'keypressfcn');
            objs = setdiff(f, objs);
            set(objs, 'keypressfcn', @(s,e)self.keypress(e));
        end

        function res = val(self, fieldname)
            % val - Helper function for processing edit boxes
            %
            %   Internal function for grabbing the value of an edit box as a
            %   number. This shortens all the getter function bodies
            %   significantly.
            %
            % USAGE:
            %   res = self.val(fieldname)
            %
            % INPUTS:
            %   fieldname:  String, Name of the field of the handles structure
            %               that we want to grab the value of.

            res = str2double(get(self.Handles.(fieldname), 'string'));
        end

        function setVal(self, fieldname, value)
            % setVal - Helper function for setting editbox values
            %
            %   Internal function for specifying the value of an input
            %   parameter that is shown in an edit box.
            %
            % USAGE:
            %   self.setVal(fieldname, value)
            %
            % INPUTS:
            %   fieldname:  String, Name of the field of the handles
            %               structure that we want to set the value of
            %
            %   value:      Number, Numeric value to set the value to

            set(self.Handles.(fieldname), 'string', num2str(value))
            notify(self, 'StateChanged');
            self.preview();
        end
    end
end
