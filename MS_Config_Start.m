%ATH IOANNIS 1444
% The following script is used to study the Diffraction Grating on a
% Metasurface for a Spherical Wave

clearvars; close all; clc;

% ----------- Users preffered General/Metasurface Parameters --------------

c0 = 3e8; % [m/s] Light speed in vacuum
wl = 0.0107; % [m] Wavelength
M =32; % # of rows
N = 32; % # of columns
Du = 0.304*wl; % [m] unit cell size
freq = c0/wl; % [Hz] Frequency Requaired for Power_Rad_final_Func

MetaSurf.ScatPat_Exp =0; % [.] exponent in expression: cos(theta)^n
MetaSurf.Rp_mnNF = [];

% 3D coordinates for the center of the wave front & Amplitude
Illumin.xyz = [ 0 , 0 ,  Inf ]; % [m] wavefront center
Illumin.Inc_Ampli = ones(M,N); % ones(M,N) = Even Amplitude / Anything else for Uneven[.]

% for the FAR-field we need the following variables to be set:
% For the Direction Of Departure we can put the sets of theta,phi angles
% i.e. if we put 1 set of theta,phi we have beam steering and from 2 sets &
% up we will have beam splitting to the given sets of directions. The sets
% must be placed in rows!
% Scatt.DOD = [  deg2rad(60) , deg2rad(0)]; % [rad] theta,phi angles
Scatt.DOD = [];
    %deg2rad(30), deg2rad(270)];
    %deg2rad(70), deg2rad(120);
    %deg2rad(50), deg2rad(260)]; % this no-scat-dir produces a default pattern

% for the NEAR-field we need the following variables to be set:
% PoF follows the same principle as the DOD matrix
Scatt.PoF = [0,0,10*wl];
                % 0 , 0 ,7*wl ]; % (Point of Focus) a 3x1 vector that decides where the Spotdiagram will be vertical from the MS
Scatt.IPD =0; % (Image Point Diagram FOR NEAR FIELD ONLY) 0 = XY Plane Spot-Diag / 1 = XZ Plane Field Plot
Scatt.HFP = 0; % 0= HFP / 1 = IHFP
% --------------------- Calculations & Initializations for Param ---------------------

MetaSurf.Du = Du; % [m] cell-width (square)
MetaSurf.NumCellsXY = [ M N ]; % [.] number of cells in each dimension
MetaSurf.Refl_Ampli = ones(M,N); % [.]
MetaSurf.Rp_mnNF = [];
MetaSurf.Ra_mnNF = [];
MetaSurf.DoPlot = 0;

% ----------- Call function to calc. scat-patt -----------
[E_wanted,A_FP, P_PF ]=HFP_ScatPat_SpotDiag_Main( freq , MetaSurf , Illumin , Scatt );
EWNRM = E_wanted / sqrt(sum(sum(abs(E_wanted).^2)));

%%
% --------------------- Gerchberg-Saxton Algorithm -------------------------
% we start the algorithm by holding onto the Amplitude profile of a
% focusing metasurface and adding a randomized phase profile to gradualy
% create the Airy disk we requaire in the XY radiation diagram
    
Illumin.Inc_Ampli = ones(M,N); % <------------- normal incidence plane wave
Illumin.xyz = [0,0,Inf]; % <--- spherical source at infinity == plane wave
MetaSurf.Rp_mnNF = rand(M,N)*2*pi; % <----- start with random phase

% >>> We're looking for the MS phase-profile for the following |Escat| <<<<
EWNRM = abs(EWNRM);

% You can load an image (gotta be 32x32 JPG)
a = sum( double( imread( 'pixil-frame-0.png' , 'png' ) ) , 3 );
a = abs(a)/max(abs(a(:)));
EWNRM = a;

% ----------------------------------------------------

% Ensure the E_wanted normalized (so that its total energy is ==1)
EWNRM = EWNRM / sqrt( sum(   abs(EWNRM(:)).^2 ) );

% Lets do some GSA iteraions
NiGSA = 30; % number of GSA iterations to do
TLR = 0.01; % tolerance, to see if GSA got us the hologram we wanted
subplot(2,5,1)
imagesc( abs(EWNRM).^2 ); title('|E_{requested}|^2'); 
axis equal tight off;

for n=1:NiGSA
           
    % Go from MS to IP (HFP). Get the E_scat there. 
    Scatt.HFP = 0; % Do HFP
    MetaSurf.Ra_mnNF = ones(N,M); % <----- MS-magn is assumed flat
    [E2, A2_prof ,Ph2_prof ]=HFP_ScatPat_SpotDiag_Main( freq , MetaSurf , Illumin , Scatt );
    E2NRM =  E2 / sqrt(sum(sum(abs(E2).^2))); 
    
    % >>>> Does the E_scat look like the hologram we wanted??
    
    % First, let's do some computations E.g., calc a weighted measure of the 
    % "error" (1-similarity). We first compute a normalized
    % difference per-pixel, and then weight that to get the average.
    myWeight = abs(EWNRM).^2; % normalized |E_wanted|^2
    compareDiff = (abs(EWNRM).^2 - abs(E2NRM).^2)./( abs(EWNRM).^2 + eps ); 
    RelErr(n) =  sum(sum(abs(compareDiff).*myWeight)) / sum(sum(myWeight)) ; 
    fprintf( ' ** Error in acquired hologram = %e\n' , RelErr(n) );
    if RelErr(n)<TLR
        disp(' !! GSA is successful! Exiting')
        figure;
        imagesc( 180/pi*(MetaSurf.Rp_mnNF+pi) ); colorbar
        colormap("jet")
        title( ' This is the MS phase-profile for the hologram' );
        return;
    elseif n>5 && std(diff(RelErr(n:-1:n-5))) < TLR/10
        disp(' ** GSA has converged (no further improvement possible)')
        figure;
        imagesc( 180/pi*(MetaSurf.Rp_mnNF+pi) ); colorbar
        colormap("jet")
        title( ' This is best MS phase-profile the GSA could get' );
        return;    
    end
       
    % Secondly, let's also plot both, for visual inspection:
    if n<=9
        spi = n+1;
    else
        spi = 10;
    end
    subplot(2,5,spi)    
    ETP = abs(E2NRM).^2;
    ETP = ETP./max(ETP(:));
    imagesc( ETP ); 
    title( sprintf( ' n_{GSA}=%d : err %4.2f%%' , n , RelErr(n)*100 ) );
    axis equal tight off;
    
    % OK, if you're here, it means that the hologram in the IP is no good.
    % So, keep the E_scat phase, replace the E_scat amplitude with what you
    % wanted (from the hologram), and IHFP back to the MS plane. Once
    % there
    
    % Go back to MS
    MetaSurf.Rp_mnNF = angle(E2); % <----- phase(E2) is now the phase(MS) for IHFP
    MetaSurf.Ra_mnNF = EWNRM; % <----- ampl(MS) is the wanted |E| for IHFP
    Scatt.HFP = 1; % 1 ---> Do IHFP
    [E3, A3_prof ,Ph3_prof ]=HFP_ScatPat_SpotDiag_Main( freq , MetaSurf , Illumin , Scatt );
    MetaSurf.Rp_mnNF = angle(E3); % <---------------- replace MS phase
    
    fprintf( ' ---------- GSA step #%d iteration end ----------\n ' , n );
    if n == NiGSA
        figure;
        imagesc( 180/pi*(MetaSurf.Rp_mnNF+pi) ); colorbar
        colormap("jet")
        title( ' This is final MS phase-profile when max N_{GSA} reached' );
    end
end
