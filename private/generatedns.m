function dns = generatedns(varargin)
    % generatedns - Generate DENSE data according to specific parameters
    %
    %   This function uses the specified parameters to generate simulated
    %   phase images as well as all other necessary information to perform
    %   DENSE analysis. The resulting structure can be stored in a .dns
    %   file and be loaded directly into DENSEanalysis software for strain
    %   calculation.
    %
    %   The input parameters are as follows:
    %
    %   InnerRadius:    Scalar, Indicates the endocardial radius in mm
    %
    %   OuterRadius:    Scalar, Indicates the epicardial radius in mm
    %
    %   RadialProfile:  [N x 2] Used to determine the amount of radial
    %                   contraction present transmurally. The first column
    %                   indicates the transmural position between 0 and 1
    %                   where 0 is the endo and 1 is the epi. The second
    %                   column indicates the maximum radial displacement at
    %                   each of these transmural locations (in mm)
    %
    %   TwistProfile:   [N x 2] Used to determine the amount of twist at
    %                   each transmural layer of the myocardium. The first
    %                   column indicates the transmural position between 0
    %                   and 1 where 0 is the endo and 1 is the epi. The
    %                   second column indicates the maximum twist at each
    %                   of these transmural locations (in radians)
    %
    %   PixelSpacing:   [1 x 2] Array, Indicates the vertical and
    %                   horizontal pixel spacing (in mm/px)
    %
    %   Frames:         Integer, Number of cardiac phases to create
    %
    %   Padding:        Scalar, Indicates the amount of padding to include
    %                   in the phase images around the epicardium
    %
    %   DENC:           Scalar, Displacement encoding value (in cyc/mm)
    %
    % USAGE:
    %   dns = generatedns(params)
    %
    % INPUTS:
    %   params: Param/Value Pairs, As specified above
    %
    % OUTPUTS:
    %   dns:    Struct, Structure similar to that contained in the .dns
    %           files from DENSEanalysis. This structure can be saved to a
    %           file for later processing

    % Parse and verify all input parameters
    ip = inputParser();
    ip.addParamValue('InnerRadius', 30, @(x)isscalar(x) && isnumeric(x));
    ip.addParamValue('OuterRadius', 40, @(x)isscalar(x) && isnumeric(x));
    ip.addParamValue('RadialProfile', [0 0; 1 0], @(x)size(x, 2) == 2);
    ip.addParamValue('TwistProfile', zeros(2), @(x)size(x, 2) == 2);
    ip.addParamValue('PixelSpacing',[1 1],@(x)isnumeric(x) && numel(x)==2);
    ip.addParamValue('Frames', 17, @(x)x >= 0 && round(x) == x);
    ip.addParamValue('Padding', 0.15, @(x)isscalar(x) && isnumeric(x));
    ip.addParamValue('DENC', 0.1, @(x)isscalar(x) && isnumeric(x));
    ip.parse(varargin{:});

    RadialProfile   = ip.Results.RadialProfile;
    TwistProfile    = ip.Results.TwistProfile;
    InnerRadius     = ip.Results.InnerRadius;
    OuterRadius     = ip.Results.OuterRadius;
    PixelSpacing    = ip.Results.PixelSpacing;
    DENC            = ip.Results.DENC;
    nFrames         = ip.Results.Frames;
    PADDING         = ip.Results.Padding;

    RADIAL_SAMPLES  = 1000;

    % +/- width
    xtra = OuterRadius * PADDING;

    % Ranges in millimeters
    xrange = ceil((OuterRadius + xtra) / PixelSpacing(2));
    xrange = -xrange:xrange;
    yrange = ceil((OuterRadius + xtra) / PixelSpacing(1));
    yrange = -yrange:yrange;

    %---------------------------------------------------------------------%
    %                Determine time evolution of Contours                 %
    %---------------------------------------------------------------------%
    epiRadiusOverTime  = interp1([1 2 3], ...
                                 OuterRadius-[0 RadialProfile(end,2) 0],...
                                 linspace(1,3,nFrames),'spline');

    endoRadiusOverTime = interp1([1 2 3],...
                                 InnerRadius-[0 RadialProfile(1,2) 0],...
                                 linspace(1,3,nFrames),'spline');

    nRadialStops    = size(RadialProfile, 1);
    nAngularStops   = size(TwistProfile, 1);

    radii   = nan(nRadialStops, nFrames);
    thetas  = nan(nAngularStops, nFrames);

    % Setup radial splines
    for stops = 1:nRadialStops
        radii(stops,:) = interp1([1 2 3],...
                                 [0; RadialProfile(stops,2); 0], ...
                                 linspace(1,3,nFrames),'spline');
    end

    for stops = 1:nAngularStops
        thetas(stops,:) = interp1([1 2 3],...
                                  [0; TwistProfile(stops,2); 0],...
                                  linspace(1,3,nFrames),'spline');
    end

    % Create endo and epicardial boundaries over the cardiac cycle
    t = linspace(0, 2*pi, 13); t(end) = []; t = t(:);

    % Anonymous function for generating contours (in pixel coordinates)
    createContour = @(r)[r .* cos(t) ./ PixelSpacing(2),...
                         r .* sin(t) ./ PixelSpacing(2)];

    endos   = arrayfun(createContour, endoRadiusOverTime, 'uni', 0);
    epis    = arrayfun(createContour, epiRadiusOverTime, 'uni', 0);

    %---------------------------------------------------------------------%
    %                Create Magnitude Images for Display                  %
    %---------------------------------------------------------------------%

    nX = numel(xrange);     midx = (nX / 2) + 0.5;
    nY = numel(yrange);     midy = (nY / 2) + 0.5;

    midxmm = midx .* PixelSpacing(2);
    midymm = midy .* PixelSpacing(1);

    % Initialize all images to zeros
    [MagnitudeImages, Xunwrap, Yunwrap] = deal(zeros([nY, nX, nFrames]));

    % Convert from pixels to millimeters
    [XX,YY] = meshgrid(xrange, yrange);
    XX = XX .* PixelSpacing(2);
    YY = YY .* PixelSpacing(1);

    % Shift endos and epis so that they are centered on the image
    func    = @(x)bsxfun(@plus,x,[midx, midy]);
    endos   = cellfun(func, endos, 'uni', 0);
    epis    = cellfun(func, epis, 'uni', 0);

    for frame = 1:nFrames

        % Generate spline for the radius and angle over the myocardium
        % Determine change in radius and angle throughout the wall
        rr = linspace(0, 1, RADIAL_SAMPLES);
        DR = spline(RadialProfile(:,1), radii(:,frame), rr);
        DT = spline(TwistProfile(:,1), thetas(:,frame), rr);

        % Parameterization of initial myocardium
        r0 = linspace(InnerRadius, OuterRadius, RADIAL_SAMPLES);

        % This is the RADIAL position of these material points across the
        % wall at the current frame
        rnew = r0 - DR;

        % Create a spline for displacements based upon the CURRENT position
        rpp = spline(rnew, DR);
        tpp = spline(rnew, DT);

        % Generate the pixel mask
        t = linspace(0, 2*pi, 13); t(end) = [];
        xy = [cos(t(:)), sin(t(:))];

        endoxy = endoRadiusOverTime(frame) .* xy;
        epixy  = epiRadiusOverTime(frame) .* xy;

        spacing = 0.5 * PixelSpacing(1);

        % This is very inefficient but that's how DENSEanalysis does it so
        % we need to ensure it is done the same way so the masks match
        seg = clinesegments(endoxy, true, true(12,1), false(12,1), spacing);
        encrv = cat(1, seg{:});
        tf = [true; all(encrv(2:end,:)~=encrv(1:end-1,:),2)];
        encrv = encrv(tf,:);

        seg = clinesegments(epixy, true, true(12,1), false(12,1), spacing);
        epcrv = cat(1, seg{:});
        tf = [true; all(epcrv(2:end,:)~=epcrv(1:end-1,:),2)];
        epcrv = epcrv(tf,:);

        [inep,onep] = inpolygon(XX, YY, epcrv(:,1), epcrv(:,2));
        [inen,onen] = inpolygon(XX, YY, encrv(:,1), encrv(:,2));

        MagnitudeImages(:,:,frame) = (inep & ~inen) | onep | onen;

        % Grab the pixel centers to see what values they have (px coords)
        [Y0, X0] = find(MagnitudeImages(:,:,frame));

        % Convert to millimeters
        X0 = X0 .* PixelSpacing(2);
        Y0 = Y0 .* PixelSpacing(1);

        % For now we only have radial gradients
        [~,R] = cart2pol(X0 - midxmm, Y0 - midymm);

        dr = ppval(rpp, R);
        dt = ppval(tpp, R);

        % dr here is positive meaning that we need to subtract to get the
        % right thing
        [TH,R] = cart2pol(X0 - midxmm, Y0 - midymm);

        % Convert dr, dt to dx, dy
        Rnew = R + dr;
        Tnew = TH + dt;

        [Xnew,Ynew] = pol2cart(Tnew, Rnew);

        % In millimeters
        dx = Xnew - X0 + midxmm;
        dy = Ynew - Y0 + midymm;

        dx = -dx;
        dy = -dy;

        xuw = Xunwrap(:,:,frame);
        xuw(logical(MagnitudeImages(:,:,frame))) = dx;
        Xunwrap(:,:,frame) = xuw;

        yuw = xuw;
        yuw(logical(MagnitudeImages(:,:,frame))) = dy;
        Yunwrap(:,:,frame) = yuw;
    end

    % DENC is in cycles / mm
    % -1 = -2 * DENC
    Xwrap = Xunwrap .* DENC; % cyc;
    Xwrap = Xwrap + 0.5;   % Add half cycle to get to zero
    Xwrap = mod(Xwrap, 1);
    Xwrap = uint16(Xwrap * 4095);

    Ywrap = Yunwrap .* DENC; % cyc;
    Ywrap = Ywrap + 0.5;   % Add half cycle to get to zero
    Ywrap = mod(Ywrap, 1);
    Ywrap = uint16(Ywrap * 4095);

    Mag = uint16(MagnitudeImages);

    % Set the upper left hand pixel to show scaling
    Mag(1,1,:) = 6;


    %% GENERATE THE DISPLACEMENT FIELD EQUATIONS %%

    nNodes = numel(endos{1}(:,1));

    % Now create the DNS file
    iscorner = {false(1, nNodes)};
    iscorner = repmat(iscorner, [nFrames, 2]);

    isclosed = num2cell(true(nFrames, 2));
    %isclosed = repmat(isclosed, [nFrames, 2]);

    iscurved = isclosed;

    dns.roi = struct(   'Name',     'Region.Of.Interest',...
                        'Type',     'SA',...
                        'UID',      dicomuid,...
                        'SeqIndex', 1,...
                        'Position', {cat(2, epis(:), endos(:))},...
                        'IsClosed', {isclosed},...
                        'IsCurved', {iscurved},...
                        'IsCorner', {iscorner});

    %% Images %%
    dns.img = {Mag; Xwrap; Ywrap};

    template = fullfile(fileparts(mfilename('fullpath')), 'template.dns');

    tmp = load(template, '-mat');
    seq = tmp.seq(1);

    empty = repmat({''}, [nFrames, 1]);

    seq.Width       = size(Mag,2);
    seq.Height      = size(Mag,1);
    seq.Filename    = empty;
    seq.FileModDate = empty;
    seq.FileSize    = empty;
    seq.MediaStorageSOPInstanceUID = empty;
    seq.InstanceCreationTime = empty;

    seq.SOPInstanceUID = empty;
    seq.AcquisitionTime = empty;
    seq.ContentTime = empty;
    seq.TriggerTime = num2cell((1:nFrames).');
    seq.NominalInterval = num2cell(1000 * ones(nFrames,1));

    seq.InstanceNumber = num2cell((1:nFrames).');

    seq.LargestPixelValue = 4095;
    seq.WindowCenter = 2048;
    seq.WindowWidth = 2048;

    seq.NumberInSequence = nFrames;

    seq = repmat(seq, [3 1]);

    seq(1).DENSEid = 'mag.overall';
    seq(1).DENSEindex = num2cell((1:nFrames).');
    seq(1).DENSEdata = struct('Number',     nFrames,...
                              'Partition',  [1 1],...
                              'Scale',      [],...
                              'EncFreq',    [],...
                              'SwapFlag',   0,...
                              'NegFlag',    [0 0 0]);

    imagecomments = empty;
    for i = 1:nFrames
        fmt = ['DENSE overall mag - Rep:0/1 Slc:0/1 Par:0/1 Phs:%d/%d ',...
                'RCswap:0 RCSflip:0/0/0'];
        imagecomments{i} = sprintf(fmt, i-1, nFrames);
    end
    seq(1).ImageComments = imagecomments;


    seq(2).DENSEid = 'pha.x';
    seq(2).DENSEindex = num2cell((1:nFrames).');
    seq(2).DENSEdata = struct('Number',     nFrames,...
                              'Partition',  [1 1],...
                              'Scale',      1,...
                              'EncFreq',    DENC,...
                              'SwapFlag',   0,...
                              'NegFlag',    [0 0 0]);

    imagecomments = empty;
    for i = 1:nFrames
        fmt = ['DENSE x-enc pha - Scale:1.000000 EncFreq:%0.2f Rep:0/1 ',...
               'Slc:0/1 Par:0/1 Phs:%d/%d RCswap:0 RCSflip:0/0/0'];
        imagecomments{i} = sprintf(fmt, DENC, i-1, nFrames);
    end
    seq(2).ImageComments = imagecomments;


    seq(3).DENSEid = 'pha.y';
    seq(3).DENSEindex = num2cell((1:nFrames).');
    seq(3).DENSEdata = struct('Number',     nFrames,...
                              'Partition',  [1 1],...
                              'Scale',      1,...
                              'EncFreq',    DENC,...
                              'SwapFlag',   0,...
                              'NegFlag',    [0 0 0]);

    imagecomments = empty;
    for i = 1:nFrames
        fmt = ['DENSE y-enc pha - Scale:1.000000 EncFreq:%0.2f Rep:0/1 ',...
               'Slc:0/1 Par:0/1 Phs:%d/%d RCswap:0 RCSflip:0/0/0'];
        imagecomments{i} = sprintf(fmt, DENC, i-1, nFrames);
    end
    seq(3).ImageComments = imagecomments;

    dns.seq = seq;

    dns.dns = struct('Name',        'Simulation.Data',...
                     'UID',         dns.roi.UID,...
                     'Type',         'xy',...
                     'MagIndex',     [1 1 NaN],...
                     'PhaIndex',     [2 3 NaN],...
                     'Number',       nFrames,...
                     'PixelSpacing', PixelSpacing,...
                     'Scale',        [1 1 NaN],...
                     'EncFreq',      [DENC DENC NaN],...
                     'SwapFlag',     0,...
                     'NegFlag',      [0 0 0]);
end
