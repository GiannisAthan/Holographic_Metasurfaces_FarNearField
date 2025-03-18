%ATH IOANNIS 1444

function [E,Amp_prof,Ph_prof] = HFP_ScatPat_SpotDiag_Main( freq , MetaSurf , Illumin , Scatt )

% FUNCTION [E,the,phi] = HFP_ScatPat_Func_Sph( freq , MetaSurf , Illumin , Scatt )
%
% ===== Inputs =====
% Check the "nargin==0" block for values of these arguments
%  * frequency [Hz] scalar
%  * MetaSurf -- structure defining the MS properties, with fields:
%     .duc [m] scalar // cell-width (square)
%     .NumCellsXY [.] 1x2 vector // # of cells in Ny/rows and Mx/columns
%     .Refl_Ampli [.]  Ny-by-Mx matrix // reflection-coeff amplitude profile
%     .Refl_Phase [rad] Ny-by-Mx matrix // reflection-coeff phase profile
%     .Cell_ScatPat_Exp [.] scalar // exponent in expression: cos(theta)^n
%  * Illumin -- structure defining the illumination/incident wave properties
%     .Type [.] Illumination/source type: 0=Plane wave, 1=Spherical
%     .DOA [rad] 1x2 vector // Direction-of-Arrival for plane wave [theta,phi]
%     .xyz [m] 1x3 vector // coordinates to the spherical wave center
%  * Scatt -- structure defining the Scattered/Reflected wave properties
%     .DOD [rad] nx2 vector // Direction-Of-Departure for plane wave [theta,phi] the n parameter represents
%  the amount of split rays we want to have in the event of Beam splitting (in the event of beam steering the vector is 1x2)
%     .scatt_Type [.] Type of scattering // 0 = Random MS Config , 1 = Beam Steering , 2 = Beam Splitting
%     .field_Dist [.] Defining Distance field for the Scattering Diagram // 0 = Far-field , 1 = Near-field
% ===== Outputs =====
%  * E -- scattering pattern (E-field complex amplitude)
%  * [the,phi] -- matrices from meshgrid, with the directions corresponding
%       to the scattering patter. They are: theta=0:SPAR:90, phi=0:SPAR:360


% -------------------------------------------------------------------------
% ============== Parameter Initialization WITHOUT Inputs ==================
% -------------------------------------------------------------------------

if nargin == 0
    clearvars;
    close all;  clc;
    disp([' Starting to Calculate and Produce the Radiation Diagram for a MxN MS ' ...
        'using default parameters ']);
    
    freq = 1e12; % [Hz] Operating frequency
    wl = 3e8 / freq;
    
    % Metasurface Standard Parameters
    MetaSurf.Du = wl/5; % [m] unit-cell size (square)
    MetaSurf.NumCellsXY = [ 40 40 ]; % [.] number of cells in each dimension
    MetaSurf.ScatPat_Exp = 0; % [.] exponent in expression: cos(theta)^n
    
    MetaSurf.Rp_mnNF = zeros(MetaSurf.NumCellsXY(1)*MetaSurf.NumCellsXY(2));
    Metasurf.Ra_mnNF = ones(MetaSurf.NumCellsXY(1)*MetaSurf.NumCellsXY(2))
    MetaSurf.DoPlot = 0;
    
    % Apartures' Phase and Amplitude
    MetaSurf.Refl_Ampli = ones(  MetaSurf.NumCellsXY(1) , MetaSurf.NumCellsXY(2) );
    
    % Illumination Angles / Angles of Arrival
    Illumin.Inc_Ampli = ones(  MetaSurf.NumCellsXY(1) , MetaSurf.NumCellsXY(2) );
    Illumin.DOA = deg2rad([ 30 , 90 ]); % [rad] theta+phi angles (plane waves)
    Illumin.Type = 0; % Source/Illumination type: 0=Plane, 1=Spherical
    
    % Scattered Wave Direction Of Departure
    Scatt.scatt_Type = 0; % 0 -> Random MS Phase Config for random Beam Steering
    Scatt.DOD=[0];
    Scatt.PoF = [];
    Scatt.HFP = 0;
end

% -------------------------------------------------------------------------
% ================== Function Parameter Preperation =======================
% -------------------------------------------------------------------------

% ------------- Metasurface Parameter Prep -------------
Du = MetaSurf.Du; % [m] unit-cell size (square)
N = MetaSurf.NumCellsXY(1); % y-dim (up/down, #-of-row)
M = MetaSurf.NumCellsXY(2); % x-dim (left/right, #-of-column)
Ra_mnFF = MetaSurf.Refl_Ampli ; % [.] amplitude profile

n = MetaSurf.ScatPat_Exp ; % exponent for lateral directions of scattered wave
DoPlot = MetaSurf.DoPlot; % Variable to not print = 0 or print = else the plots

% Converting the xyz coordinates to spherical
ro = sqrt(Illumin.xyz(1)^2 + Illumin.xyz(2)^2 + Illumin.xyz(3)^2);
phi_i = acos(Illumin.xyz(3)/ro);
theta_i = sign(Illumin.xyz(2)) * acos(Illumin.xyz(1) / ( sqrt(Illumin.xyz(1)^2 + Illumin.xyz(2)^2) ) );

% We use the following if to check the distance of the wave center in
% relevance to the metasurface. If the distance is greater than 10^5 times
% the wavelength of the wave we consider it a plane wave, else its spherical.
if isinf(ro)
    
    disp ('Incident Wave is in Far-Field, as a Plane Wave')
    Illumin.Type = 0; % Source/Illumination type: 0=Plane, 1=Spherical
    
    Illumin.DOA = [ 0 , 0 ]; % (MUST BE CALLED FOR PLANE WAVE) [rad] theta+phi angles
else
    
    disp ('Incident Wave is in Near-Field, as a Spherical Wave')
    Illumin.Type = 1; % Source/Illumination type: 0=Plane, 1=Spherical
    
    if Illumin.xyz(1) ==0 && Illumin.xyz(2) == 0
        Illumin.DOA = [ 0 , 0 ]; % (MUST BE CALLED FOR PLANE WAVE) [rad] theta+phi angles
    end
end

Scatt.scatt_Type = size(Scatt.DOD,1); % [.] 1 = Beam Steering / 2 >= Beam Splitting

% ------------- Illumination Parameter Prep -------------
ill_Type = Illumin.Type; % Source/Illumination type: 0=Plane, 1=Spherical

if ill_Type == 0
    theta_i = Illumin.DOA(1); % [rad] theta incident/illum (for plane)
    phi_i = Illumin.DOA(2); % [rad] phi incident/illum (for plane)
elseif ill_Type ==1
    Sph_Cntr = Illumin.xyz; % [m] sphere center (for spherical)
else
    error(['The value of the ill_Type varaible must be either 0 for PW' ...
        'or 1 for Spherical Wave.'])
end

%-------------- Scattering/ Reflection Parameters ----------------------
% Defining the Distance of the plotting Diagram in refrence to the MS
Scatt_Type = Scatt.scatt_Type; % 0 = Random MS Config , 1 = Beam Steering , 2 = Beam Splitting

% -------------------------------------------------------------------------
% =================== Basic Parameter Initialization ======================
% -------------------------------------------------------------------------
wl = 3e8 / freq;
k0 = 2*pi/wl; % [rad/m] free-space wavenumber

%creating MS grid and Parameters
[mm,nn] = meshgrid(1:M,1:N); % index-array in MS rows/columns
dxk = k0 * Du; % unit cell size, x-dimension*kappa (d*k0=2*pi * d/wl)
dyk = k0 * Du; % unit cell size, y-dimension*kappa

% Ms grid parameters used for the Spherical Wave & Near-Field Calculations
[xuc,yuc] = meshgrid( (-(M-1)/2:(M-1)/2)*dxk , (-(N-1)/2:(N-1)/2)*dyk ) ; % [.]
zuc = k0*zeros(size(xuc)); % [.] N-by-M size (rs has M*N cells)

% -------------------------------------------------------------------------
% ================ Incident Wave/Source Parameters ========================
% -------------------------------------------------------------------------

% For simplification purposes we use a Plane Wave for the Incident Wave and
% we set the amplitude to 1.

% Source/Illumination type: 0=Plane, 1=Spherical
if ill_Type == 0
    % Direction of Source/Incident Wave , for each unit-cell
    % For Plane-Wave (PW) incidence, theta & phi are constant across the MS
    theSmn = zeros(M,N) + theta_i; % (rad)
    phiSmn = zeros(M,N) + phi_i; % (rad)
    
    % Amplitude and Phase of Incident Plane Wave
    Ia_mn = ones(M,N) ;
    Ip_mn = ( dxk*(mm-1).*cos(phiSmn) - dyk*(nn-1).*sin(phiSmn) ).*sin(theSmn);
    IncDirec = rad2deg([theta_i, phi_i]); % [deg] Calling theta,phi DOA
    
elseif ill_Type == 1
    
    % If the Spherical Waves center is in the far-field the Wave reaches
    % the MS as a Plane Wave so we terminate the Function.
    if Sph_Cntr(3) < 3*wl
        error('Spherical Wave Center too close to m/s. Stopping function execution.')
    end
    
    xyzks = k0*Sph_Cntr; % [.] spher-source xyz, normed to kappa0
    
    % For Spherical-Wave (SW, point source), theta & phi are rigorously
    % calculated from geometric formulas.Using the meshgrid 104-106
    
    % 1. Transform source coords
    xsrc = xyzks(1);
    ysrc = xyzks(2);
    zsrc = xyzks(3);
    % 2. cartesians dists of source from each unit-cell center
    xd_sfuc = xsrc - xuc;
    yd_sfuc = ysrc - yuc;
    zd_sfuc = zsrc - zuc;
    % 3. transform to sphericals
    dist_suc = sqrt(xd_sfuc.^2 + yd_sfuc.^2 + zd_sfuc.^2);
    
    theSmn = atan(sqrt(xd_sfuc.^2 + yd_sfuc.^2)./zd_sfuc); % (rad)
    phiSmn = atan2( yd_sfuc, xd_sfuc ); % (rad)
    
    % Point Source :
    if Illumin.Inc_Ampli == ones(M,N)
        % for even Amplitude Illumination
        disp(' Using the default EVEN Amplitude for Spherical Wave = 1 ')
        Ia_mn = ones(M,N);
    else
        % for uneven Amplitude Illumination
        disp(' Using the UNEVEN Amplitude for Spherical Wave ')
        R_cs = norm(dist_suc);
        
        % Giving E0 this value so the center of the MS has an Amplitude of 1
        E0 = 4*pi*mean(R_cs);
        
        Ia_mn = E0./40*(4*pi*dist_suc);
    end
    Ip_mn = dist_suc; % Because phase=beta*dz-->k0*distance
    IncDirec = rad2deg([theSmn, phiSmn]); % [deg] Calling theta,phi DOA for Sph;
end

% -------------------------------------------------------------------------
% ======= MS Grating Calculations / Scattered wave Parameter Prep =========
% -------------------------------------------------------------------------
% Given input phase at each unit cell (from Dir-of-Arrival)
Phase_inc = Ip_mn ;

% Depending on the Scatt_Type we have:
% 0 -> Random MS Phase Config for random Beam Steering
% 1 -> Beam Steering to the theta and phi angles that we used to call Func
% 2 -> Beam Splitting in 2 or more directions rows in DOD decide the split
% lobes
if ((Scatt_Type == 0) && isempty(Scatt.DOD) && isempty(Scatt.PoF)) || (nargin == 0)
    
    % ---------------------------------------------------------------------
    % ====== Apartures' Phase Initialization WITHOUT parameters given =====
    % ---------------------------------------------------------------------
    disp(' Creating phase profile with random noise (diffused scattering) ')
    Rp_mnFF = rand(size(nn))*2*pi;
    
elseif Scatt_Type >= 1
    
    if Scatt_Type == 1
        disp(' Performing Beam Steering ')
    elseif Scatt_Type >= 2
        disp([' Performing Beam Splitting in ', num2str(Scatt_Type), ' lobes '])
    end
    
    % Initialize Rp_mn to zeros
    Rp_mnFF = zeros(M, N);
    
    for lobe = 1:Scatt_Type
        % Get the direction of departure (DOD) for the current lobe
        theta_sL = Scatt.DOD(lobe, 1);
        phi_sL = Scatt.DOD(lobe, 2);
        
        % Calculate the output phase at each unit cell (for DOD)
        Phase_scaL = (-dxk*(mm-1).*cos(phi_sL) + dyk*(nn-1).*sin(phi_sL)).*sin(theta_sL);
        
        % Calculate the reflection phase coefficient of the MS for the current lobe
        Refl_Coef_phL = Phase_scaL - Phase_inc;
        
        % Ensure that the phase of each cell is inside the [0, 2*pi) rad
        Phi_mnL = mod(Refl_Coef_phL, 2*pi);
        
        % Accumulate the reflection phase coefficients for each lobe
        Rp_mnFF = Rp_mnFF + exp(1i*Phi_mnL);
    end
    % Take the angle of the accumulated complex reflection coefficients
    Rp_mnFF = angle(Rp_mnFF);
    Rp_mnFF = mod( Rp_mnFF , 2*pi );
    
end

% -------------------------------------------------------------------------
% === Huygens-Fresnel Principle for Radiation Diagram Far OR Near Field ===
% -------------------------------------------------------------------------
% ---------------- Calculating Far-field Scattering  ----------------------
if (~isempty(Scatt.DOD)) || ((Scatt_Type == 0) && isempty(Scatt.DOD) && isempty(Scatt.PoF))
    
    disp ([' Calculation of the complex E-field for Far-Field scattering using' ...
        ' the Huygens-Fresnel Principle for the given values '])
    
    % Angle-resolution of Scatt-Patt:
    dsp = 3; % degree-step-phi
    dst = 3; % degree-step-theta
    angspan = 0; % [deg] ang-span on lobe-focus-direction
    thef = 0; % [deg] lobe-focus theta
    phif = 0; % [deg] lobe-focus phi
    
    % Directions on "hemisphere" space. If angspace>0, it only calculates
    % around a specific direction (e.g. around the approx reflection-lobe)
    % to cut-down on simulation times.
    if angspan > 0
        phi = deg2rad( (-angspan/2:dsp:angspan/2) + phif );
        the = deg2rad( (-angspan/2:dst:angspan/2) + thef );
        
        if any( the < 0 )
            in = the<0;
            phi( in ) = phi(in)+180;
            the = abs(the);
        end
        
    else
        phi = (0:dsp:360)*pi/180 ;
        the = (0:dst:90)*pi/180 ;
    end
    
    % Directions in the "hemisphere", gridded
    [phi,the]=meshgrid(phi,the);
    
    % Scat-Patt formation:
    EFF = 0*the; % this will hold the pattern (complex)
    
    for m1 = 1:M
        for n1 = 1:N
            
            % Variable that denotes the scattering pattern of the mn-th unit
            % cell , the choice n=1 describes real-world dipolar scatterers.
            % When n==0, scatterers scatter isotropically (unrealistic).
            fmni = cos(theSmn(n1,m1))^n; % how each cell "gathers" from inc dirs
            fmns = cos(the).^n; % how each cell "diffuses" in scat dirs
            
            % complex-valued reflected/scattered profile
            R_w = Ra_mnFF(n1,m1).* fmns .* exp( 1i * Rp_mnFF(n1,m1) );
            
            % complex-valued illumination/incidence profile
            Il_w = Ia_mn(n1,m1) .* fmni .* exp( 1i * Ip_mn(n1,m1) );
            
            Il_R_w =  Il_w .* R_w ;
            
            % Array parameter
            zeta_mn = exp( +1j*dxk*m1.*sin(the).*cos(phi))...
                .*exp( -1j*dyk*n1.*sin(the).*sin(phi));
            
            EFF = EFF + Il_R_w .* zeta_mn;
        end
    end
    
    Rad_LinearFF = abs(EFF).^2/max(abs(EFF(:)).^2);
    Rad_dBFF= 10 * log10(Rad_LinearFF);
end

% ----------------- Calculating near field HFP --------------------
if ~isnan(Scatt.PoF)
    disp([' Calculation of the complex E-field regarding the Near-Field' ...
        ' Scattering using the Fresnel-Kirchhoff Principle '])
    
    % Recall the xyz coordinates of the (centers of all the) unit cells.
    % These points are the "rs".Using the meshgrid 104-106. denormalize them
    xuc = xuc / k0; % [m]
    yuc = yuc / k0; % [m]
    zuc = zuc / k0; % [m]
    
    % Assign the first Point of Focus
    pofx = Scatt.PoF(1);
    pofy = Scatt.PoF(2);
    pofz = Scatt.PoF(3);
    df =sqrt(pofx^2 + pofy^2 + pofz^2);
    
    % Check the z coordinate condition for each Point of Focus
    if df < 3 * wl
        error(['The z coordinate for the Point of Focus must be greater ' ...
            'than 3 times the wavelength of the incoming wave!'])
    end
    
    % Now defining the Image Plane in which the radiation Diagram will be
    % portrayed.
    K = M; % this is the pixel-number in one dimension (rp has K^2 "pixels")
    x1d = linspace(-1/2,+1/2, K) * M * Du; % [m]
    switch Scatt.IPD
        case 0 % xy-plane at some z ---> you should see a spot!
            % [u,v] = meshgrid(x1d); % [m] these two are K-by-K matrix
            u = xuc;
            v = yuc;
            w = pofz + 0*u; % [m] z-coord of points on image-plane, K-by-K
        case 1 % xz-plane at y=0 ---> you should see E-field diffraction
            z1d = linspace(wl, 1.5*df, K+30); % [m]
            [u,w] = meshgrid(x1d,z1d); % [m] these two are K-by-K matrix
            v = 0*u; % [m] y-coord of points on image-plane, K-by-K
        otherwise
            error(['Value Scatt.IPD must be 0 = XY Plane Spot-Diag' ...
                ' / 1 = XZ Plane Field Plot'])
    end
    
    if isempty(MetaSurf.Rp_mnNF)
        % Calculate the distance between each rs (unit cell) to the focus.
        % Its the distance between each rs on the MS (z=0).
        ro = sqrt((pofx - xuc).^2 + (pofy - yuc).^2 ); % MxN matrix
        Rp_mnNF = k0* (df - sqrt(df^2 +ro.^2));
        % Required parabolic phase to have on-axis focus spot
        Rp_mnNF = mod(Rp_mnNF, 2*pi);
    else
        Rp_mnNF = MetaSurf.Rp_mnNF;
    end
    
    if isempty(MetaSurf.Ra_mnNF)
        % Creating the aparture so the Airy Disk can be created
        % Determine the center of the unit cell array
        center_x = (M-1) / 2;
        center_y = (N-1) / 2;
        
        % Define the radius of the circle (in units of cell indices)
        radius = (M-1) / 2;
        % Calculate the distance of each unit cell from the center
        distance_from_center = sqrt((mm-1 - center_x ).^2 + (nn-1 - center_y ).^2);
        
        Ra_mnNF = ones(size(xuc)); % unitary amplitude (full reflection)
        % Set Rp_mn to zero outside the defined circle
        Ra_mnNF(distance_from_center >= radius) = 0;
    else
        Ra_mnNF = MetaSurf.Ra_mnNF;
    end
    
    % Define the vertical vector
    u_vert = [0, 0, 1];
    
    % Initialize the combined E-field
    ENF = zeros(size(u));
    
    for i = 1:M*N
        
        rs = [xuc(i), yuc(i), zuc(i)]; % rs = xyz of this cell of the MS
        
        % 3D-vectors connecting rs to all rp points
        dspx = u - rs(1); % K-by-K
        dspy = v - rs(2);
        dspz = w - rs(3);
        
        % dot_dsp_uvert = dspx*u_vert(1) + dspy*u_vert(2) + dspz*u_vert(3); % K-by-K
        norm_dsp = sqrt( dspx.^2 + dspy.^2 + dspz.^2 ); % K-by-K
        % theta_scat = acos( dot_dsp_uvert ./ ( norm_dsp * norm(u_vert) ) );
        
        % Variable that denotes the scattering pattern of the mn-th unit
        % cell , the choice n=1 describes real-world dipolar scatterers.
        % When n==0, scatterers scatter isotropically (unrealistic).
        % fmni = cos(theSmn(i)).^n; % how each cell "gathers" from inc dirs
        % fmns = cos(theta_scat).^n; % how each cell "diffuses" in scat dirs
        fmni = 1; %cos(theSmn(i)).^n; % how each cell "gathers" from inc dirs
        fmns = 1; %cos(theta_scat).^n; % how each cell "diffuses" in scat dirs
        
        % Calculating the obliquity factor
        % Psi = 0.5 * (cos(theSmn(i)) + cos(theta_scat)); % MxN matrix
        Psi = 1;
        
        % Forward or Inverse HFP near field
        if Scatt.HFP == 0
            % Forward HFP
            % Calculating the Distance related Var
            Dst_Var = exp(1i*k0.*norm_dsp)./norm_dsp; % MxN matrix
        elseif Scatt.HFP ==1
            % Inverse HFP
            % Calculating the Distance related Var
            Dst_Var = exp(-1i*k0.*norm_dsp)./norm_dsp; % MxN matrix
        end
        % complex-valued reflected/scattered profile
        R_w = Ra_mnNF(i) .* fmns .* exp(1i * Rp_mnNF(i));
        % complex-valued illumination/incidence profile
        Il_w =  Ia_mn(i).* fmni .* exp(1i * Ip_mn(i));
        Il_R_w =  Il_w .* R_w ;
        
        CF = (1/(1i*wl)) .* Il_R_w .* Psi .* Dst_Var;
        ENF = ENF + CF;% Accumulate CF into E
    end
    
    % Calculate the dB values for the plot
    Rad_LinearNF = abs(ENF).^2/max(abs(ENF(:)).^2);
    Rad_dBNF= 10 * log10(Rad_LinearNF);
    
    % Find the x-coordinate where E reaches its max value
    [max_ENF, max_idx] = max(Rad_LinearNF(:));
    [max_row, max_col] = ind2sub(size(Rad_LinearNF), max_idx);
    x_max_E = x1d(max_col);
    disp(['The x-coordinate where E reaches its max value: ', num2str(x_max_E)])
    
    % Find the first x-coordinate where E reaches its min value
    [min_ENF, min_idx] = min(Rad_LinearNF(:));
    [min_row, min_col] = ind2sub(size(Rad_LinearNF), min_idx);
    x_min_E = x1d(min_col);
    disp(['The first x-coordinate where E reaches its min value: ', num2str(x_min_E)])
end


% Calculations completed: Assign output arguments
if ~isempty(Scatt.DOD) && isempty(Scatt.PoF)
    E = EFF;
    Amp_prof = Ra_mnFF;
    Ph_prof = Rp_mnFF;
elseif isempty(Scatt.DOD) && ~isempty(Scatt.PoF)
    E = ENF;
    Amp_prof = Ra_mnNF;
    Ph_prof = Rp_mnNF;
end

if DoPlot == 0, 
    return; 
else

% -------------------------------------------------------------------------
% == Inc Plot / MS PHASE & AMPLITUDE CONF / Creating the Power Radiation ==
% -------------------------------------------------------------------------

% ----------------- Incident Uneven Amplitude config  -------------------
if Illumin.Inc_Ampli ~= ones(M,N)
    % Plot the incident amplitude distribution
    figure;
    imagesc(Ia_mn);set(gca,'YDir','Normal')
    axis equal tight;
    xlabel('("M") #cells'); ylabel('("N") #cells');
    title('Uneven Amplitude Distribution for Spherical Wave Illumination');
    set(gca,'XMinorTick','on','YMinorTick','on')
    colormap(jet);
    hcb1 = colorbar;
    caxis([0 1])
end

% -------------------------------------------------------------------------
% ----------------- Far-Field MS Config & Power Radiation Diagram ---------
% -------------------------------------------------------------------------
if (~isempty(Scatt.DOD)) || ((Scatt_Type == 0) && isempty(Scatt.DOD) && isempty(Scatt.PoF))
    
    figure('Name', 'MS Phase and Amplitude Config FF');
    
    % Phase-profile across the MS
    subplot(1,2,1);
    imagesc(fliplr(rad2deg(Rp_mnFF)));
    set(gca, 'YDir', 'Normal');
    axis equal tight;
    xlabel('("M") #cells'); ylabel('("N") #cells');
    title('\angle\Phi_{nm} (deg)');
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    colormap(gca, hsv); % Set colormap for this subplot
    caxis([0 360]);
    colormap(jet)
    hcb1 = colorbar;
    set(hcb1, 'YTick', 0:45:360);
    
    % Uneven amplitude plot across the MS
    subplot(1,2,2);
    imagesc(fliplr(abs(Ra_mnFF)));
    axis equal tight;
    set(gca, 'YDir', 'Normal');
    colorbar;
    colormap(gca, hot); % Set colormap for this subplot
    xlabel('("M") #cells'); ylabel('("N") #cells');
    title('|r_{nm}|');
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    caxis([0 1]);
    colormap(jet)
    
    
    figure('Name','Far-Field Power Radiation Diagram');
    
    imagesc(0:360, 0:90,Rad_dBFF) ; set(gca,'YDir','Normal')
    set(gca,'XTick',0:45:360) % horz-axis tick-labels
    colormap(hot);
    colorbar;
    title('Magnitude of Radiation Pattern (dB)');
    xlabel('\phi');
    ylabel('\theta');
    caxis([-30 0]); % don't show too negative values. -30 is enough for "zero"
    
    figure('Name','3D Pattern for Far-Field');
    plot_3D_Pattern_altered( rad2deg(the), rad2deg(phi), Rad_LinearFF , [3 6] ,IncDirec )
    hold on;
    
    % ------- Beam Steering/Splitting plotting the centers of each lobe --------
    if Scatt_Type >= 1
        
        subplot (1,2,2)
        % Plotting the center of the split rays
        for n1 = 1:Scatt_Type
            center_cos = rad2deg(Scatt.DOD(n1,1).*cos(Scatt.DOD(n1,2)));
            center_sin = rad2deg(Scatt.DOD(n1,1).*sin(Scatt.DOD(n1,2)));
            % Plotting white dot for the incident wave
            plot3(center_cos,center_sin,0,'*','Color','k','MarkerFaceColor','k','MarkerSize',10)
        end
    end
end
% -------------------------------------------------------------------------
% --------------------- Near-Field MS Conf & Near-Field Plots -------------
% -------------------------------------------------------------------------

if ~isempty(Scatt.PoF)
    
    figure('Name', 'MS Phase and Amplitude Config NF');
    
    % Phase-profile across the MS
    subplot(1,2,1);
    imagesc(fliplr(rad2deg(Rp_mnNF)));
    set(gca, 'YDir', 'Normal');
    axis equal tight;
    xlabel('("M") #cells'); ylabel('("N") #cells');
    title('\angle\Phi_{nm} (deg)');
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    colormap(gca, hsv); % Set colormap for this subplot
    caxis([0 360]);
    colormap(jet)
    hcb1 = colorbar;
    set(hcb1, 'YTick', 0:45:360);
    
    % Uneven amplitude plot across the MS
    subplot(1,2,2);
    imagesc(fliplr(abs(Ra_mnNF)));
    axis equal tight;
    set(gca, 'YDir', 'Normal');
    colorbar;
    colormap(gca, hot); % Set colormap for this subplot
    xlabel('("M") #cells'); ylabel('("N") #cells');
    title('|r_{nm}|');
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
    caxis([0 1]);
    
    
    if Scatt.IPD == 0
        % Plotting the spot diagram. To have a succesfull focusing point above
        % the metasurface we need a Airy disc to form on the xy plane.
        
        figure('Name', 'Near-Field Radiation Pattern in the XY Plane');
        
        subplot(1,2,1)
        imagesc( x1d/wl, x1d/wl, Rad_LinearNF ); set(gca,'YDir','normal')
        title( 'Linear |E|^2' );
        colorbar
        xlabel( 'x_{rp} / \lambda' );
        ylabel( 'y_{rp} / \lambda' );
        caxis([0 1])
        axis equal tight
        
        subplot(1,2,2)
        imagesc( x1d/wl, x1d/wl, Rad_dBNF ); set(gca,'YDir','normal')
        title( 'dB |E|^2' );
        colorbar;
        xlabel( 'x/\lambda' );
        ylabel( 'y/\lambda' );
        caxis([-30 0]); % don't show too negative values. -30 is enough for "zero"
        axis equal tight
        
        % --------------------- Near-Field Field Plot -----------------------------
    else
        % Assuming Unorm_dB and x1d are already defined in your script.
        % We will plot Unorm_dB against x1d with a constant z (df).
        
        figure('Name', 'Near-Field Radiation Pattern in the XZ Plane');
        imagesc(x1d/wl, z1d/wl , Rad_LinearNF);
        colormap(jet);
        colorbar;
        title('Linear |E|^2');
        xlabel('x/\lambda');
        ylabel('z/\lambda');
        axis xy; % Ensures the origin is at the bottom-left
        caxis([0 1]);
        
        % creating a line on the field plot to mark the z = df plane
        hold on
        for kdf = 1:length(df)
            line(xlim, df(kdf)*[1 1] / wl, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
            x_limits = xlim; % Get the current x-axis limits
            text(x_limits(2) * 0.95, df(kdf) / wl, ...
                sprintf( 'z = df_%d' ,kdf ) ,...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontSize', 10, 'Color', 'k');
        end
    end
    
end
end
end

