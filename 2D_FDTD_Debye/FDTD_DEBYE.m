clc 
close all
clear all

% Electromagnetic Constants
c = 299792458;
Eo = 8.854187817620e-12;
Uo = 1.2566370614e-6;

% FDTD Constants
% Number of Voxels in the X and Y axis
Nx = 256;
Ny = 256;

% Size of the Square
dx = 1.1e-3; % 1mm in the X-Direction
dy = 1.1e-3; % 1mm in the Y-Direction

% Time Step Between each iteration
dT = dx / ( c * sqrt(3) );

% Number of Iterations
Nt =5000;

% Define the Material
Eps = ones(Nx,Ny);
Mu = ones(Nx,Ny); % Initilise to Freespace Ur = 1

ECond = 0.*ones(Nx,Ny);
MCond = zeros(Nx,Ny);
cmpx = sqrt(-1);

Nt = 3000;
t0 = 180e-12;
spread = 60e-12;
for n = 1:Nt
	time(n) = n * dT;
	pulse(n) = exp(-0.5 * ((time(n) - t0) / spread)^2);
end
 
figure
subplot(2,1,1)
plot(time./1e6, pulse) 
ylabel('Amplitude');
xlabel('Time [uSeconds]')
title('Time-Domain Excitation Signal')
subplot(2,1,2)
[f PULSE]  = fft1(pulse, 1 / dT);
PULSE = abs(PULSE);
PULSE = PULSE / max(PULSE);
PULSE = 20.*log10(PULSE);
plot(f ./1e9, PULSE) 
title('Frequency-Domain Excitation Signal')
ylabel('Amplitude');
xlabel('Frequency [GHz]')

% Initiliase the Coefficeints
ca = zeros(Nx,Ny);
cb = zeros(Nx,Ny);
da = zeros(Nx,Ny);
db = zeros(Nx,Ny);

Y1 = zeros(Nx,Ny);
Y2 = zeros(Nx,Ny);
Y3 = zeros(Nx,Ny);
Y4 = zeros(Nx,Ny);
Y5 = zeros(Nx,Ny);
Y6 = zeros(Nx,Ny);

N1 = zeros(Nx,Ny);
N2 = zeros(Nx,Ny);
N3 = zeros(Nx,Ny);
N4 = zeros(Nx,Ny);
N5 = zeros(Nx,Ny);
N6 = zeros(Nx,Ny);
N7 = zeros(Nx,Ny);
N8 = zeros(Nx,Ny);
N9 = zeros(Nx,Ny);
N10 = zeros(Nx,Ny);
N11 = zeros(Nx,Ny);

%Ez Mode ---------
Ez = zeros(Nx,Ny);
Dz = zeros(Nx,Ny);
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);
%-----------------

% Previous Time Field Values Required for the Mur ABC
% The Current Field is at time n+1 using convention from FDTD book

% These are n  Fields
Ez_Old = zeros(Nx,Ny);
Dz_Old = zeros(Nx,Ny);

% These are n-1 Fields
Ez_Old_Old = zeros(Nx,Ny);
Dz_Old_Old = zeros(Nx,Ny);


% These are n-2 Fields
Ez_Old_Old_Old = zeros(Nx,Ny);
Dz_Old_Old_Old = zeros(Nx,Ny);

% Debye Parameters
Eps_static = ones(Nx,Ny);
Eps_m      = 1 .* ones(Nx,Ny);
Eps_inf    = 1. * ones(Nx,Ny);
Tau1       = 1e-16 .* ones(Nx,Ny);
Tau2       = 1e-15 .* ones(Nx,Ny);
ECond      = zeros(Nx,Ny);

fprintf('Creating Tissue Map\n');
fflush(stdout);

load Data.mat
Materials=[
1 105.174 1.0008 56.2477 1.00E-011 1.28E-009 0.4936 103.9682 % Grey Matter
2 62.5297 7.5771 33.5359 1.28E-011 9.38E-010 0.2977 35.7652  % White Matter
3 77.3415 7.7648 30.3619 1.08E-011 1.20E-009 1.5376 20.6899  % Blood
4 84.9341 1.3405 18.4736 6.87E-012 1.68E-009 2.2061 8.995    % CSF
5 16.0833 4.0694 8.02240 1.71E-011 4.61E-010 0.0642 0.85788  % Skull
6 61.5213 21.776 35.8045 2.09E-011 6.18E-010 0.6183 6.928    % Dura
7 6.90520 1.1793 3.23190 1.81E-011 1.13E-009 0.0261 0.020176 % Fat
8 62.7337 20.049 43.4722 1.96E-011 1.10E-009 0.3849 6.7656]; % Skin

Model = zeros(Nx,Ny);

ind = 1;
for i=1:256
        for j=1:256
            
            mat = Slice(i,j);
            index = mat;

            if (index == 98)
                Model(i,j) = 7; % Fat Tissue
            end
            if (index  == 4)
                Model(i,j) = 5; % Skull
            end
            if (index ==2 )
                Model(i,j) = 4; % Spinal Fluid
            end
            if (index == 117) || (index == 89)|| (index ==124) || (index == 95)
                Model(i,j) = 1; % Grey Matter
            end
            if (index == 83)
                Model(i,j) = 2; % White MAtter
            end
            
            if (index == 84)
                Model(i,j) = 6; % Dura Layer
            end
            
            if (index == 23)
                Model(i,j) = 3; % Blood
            end
        end
end


for i = 1 : Nx
    for j = 1 : Ny
	index = find(Model(i,j) == Materials(:,1));
	if (index ~= 0) 
     		Eps_inf(i,j)    = Materials(index,3);
     		Eps_static(i,j) = Materials(index,2);
     		Eps_m(i,j)      = Materials(index,4);
     		Tau1(i,j)       = Materials(index,5);
     		Tau2(i,j)       = Materials(index,6);
     		ECond(i,j)      = Materials(index,7);
    	end
    end
end


Eps = Eps .* Eo;
Mu =  Mu  .* Uo;

fprintf('Building Debye Constants\n');
fflush(stdout);

for i = 1 : Nx
    for j = 1 : Ny
        Y1(i,j) = Eo * Eps_static(i,j) + ECond(i,j) * (Tau1(i,j) + Tau2(i,j));
        Y2(i,j) = Eo * Eps_inf(i,j) * Tau2(i,j) + Eps_m(i,j) * Eo * Tau1(i,j) - Eps_m(i,j) * Eo * Tau2(i,j) + Eo * Eps_static(i,j) * Tau2(i,j) + ECond(i,j) * Tau1(i,j) * Tau2(i,j);
        Y3(i,j) = Tau1(i,j) * Tau2(i,j) * Eo * Eps_inf(i,j);
        Y4(i,j) = (Tau1(i,j) + Tau2(i,j)) ;
        Y5(i,j) = 1 * Tau1(i,j) * Tau2(i,j)  ;
        Y6(i,j) = ( ECond(i,j)/2 + Y1(i,j)/(3*dT) + Y2(i,j) /(2 * dT * dT) + Y3(i,j) / (dT*dT*dT) ) ;
 
        N1(i,j) = 3*dT*Y6(i,j);
        N2(i,j) = Y4(i,j);
        N3(i,j) = 2 * dT^2 * Y6(i,j);
        N4(i,j) = Y5(i,j);
        N5(i,j) = dT^3*Y6(i,j);
        N6(i,j) = (ECond(i,j)/2 - Y1(i,j)/(3*dT) + Y2(i,j)/(2*dT^2) - Y3(i,j)/(dT^3));
        N7(i,j) = Y6(i,j);
        N8(i,j) = Y2(i,j)/(2*dT^2) + 3 *Y3(i,j)/(dT^3);
        N9(i,j) = Y6(i,j);
        N10(i,j) =  Y2(i,j)/(2*dT^2) -3*Y3(i,j)/(dT^3);
        N11(i,j) = Y6(i,j);
        
        ca(i,j) = ( 2 * Eps(i,j) -  dT * ECond(i,j) ) / ( 2 * Eps(i,j) + dT * ECond(i,j) );
        cb(i,j) = ( 2 * dT ) / ( (2 * Eps(i,j) + dT * ECond(i,j)) * dx );
        
        da(i,j) = ( 2 * Mu(i,j) -  dT * MCond(i,j) ) / ( 2 * Mu(i,j) + dT * MCond(i,j) );
        db(i,j) = ( 2 * dT ) / ( (2 * Mu(i,j) + dT * MCond(i,j)) * dx );
    end
end

% We No longer need this 2-D arrays. 
clear Y1 Y2 Y3 Y4 Y5

% Source Position
Sx = 57;
Sy = 128;
%Ez Mode ---------
Ez = zeros(Nx,Ny);
Dz = zeros(Nx,Ny);
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);
%-----------------


% Previous Time Field Values Required for the Mur ABC
% The Current Field is at time n+1 using convention from FDTD book

% These are n  Fields
Ez_Old = zeros(Nx,Ny);
Dz_Old = zeros(Nx,Ny);

% These are n-1 Fields
Ez_Old_Old = zeros(Nx,Ny);
Dz_Old_Old = zeros(Nx,Ny);


% These are n-2 Fields
Ez_Old_Old_Old = zeros(Nx,Ny);
Dz_Old_Old_Old = zeros(Nx,Ny);

Freq = 1e9;
Ez_cmpx = zeros(Nx, Ny);

fprintf('Running Simulation\n');
fflush(stdout);

% Main Time Step
for n = 1 : Nt
    msg = sprintf('Computing Time Step: %i', n); 
    disp(msg)
     fflush(stdout);
    i = 1 : Nx - 1;
    j = 1 : Ny - 1;

    Hx(i,j) = da(i,j) .* Hx(i,j) - db(i,j) .*  ( Ez(i,j+1) - Ez(i,j) );
    Hy(i,j) = da(i,j) .* Hy(i,j) + db(i,j) .*  ( Ez(i+1,j) - Ez(i,j) );

    i = 2 : Nx;
    j = 2 : Ny;
    Dz(i,j) = Dz(i,j) + dT * ( (Hy(i,j) - Hy(i-1,j))/dx -  (Hx(i,j) - Hx(i,j-1))/dy);
    Dz(Sx,Sy) = Eo * pulse(n);
   
    i = 1:Nx;
    j = 1:Ny;
    Ez(i,j) = ((Dz(i,j) - Dz_Old_Old_Old(i,j) ) ./ N1(i,j) ) + ( N2(i,j) .* (Dz(i,j) - Dz_Old(i,j) - Dz_Old_Old(i,j) + Dz_Old_Old_Old(i,j) ) ./ N3(i,j)) + ( N4(i,j) .* ( Dz(i,j) - 3 .* Dz_Old(i,j) + 3  .* Dz_Old_Old(i,j) - Dz_Old_Old_Old(i,j) )./N5(i,j) ) + (-N6(i,j) .* Ez_Old_Old_Old(i,j) ./ N7(i,j)) + (N8(i,j) .* Ez(i,j) ./ N9(i,j)) + ( N10(i,j) .* Ez_Old_Old(i,j) ./ N11(i,j));
     
	% Perform DFT
        Ez_cmpx(:,:) = Ez_cmpx(:,:) +  exp(-cmpx * 2 * pi * Freq * n * dT) .* Ez;
   
    % E-FIELD ABC --------------------------------------------------------- 
   % Mur's ABC Boundary
   Mur = (c * dT - dx) / (c * dT + dx);
    
    for j = 2 : Ny - 1
        MurCa = Ez_Old(1,j+1) - 2 * Ez_Old(1,j) + Ez_Old(1,j-1);
        MurCb = Ez_Old(2,j+1) - 2 * Ez_Old(2,j) + Ez_Old(2,j-1);
        EQ1 = Mur * ( Ez(2,j) + Ez_Old_Old(1,j) );
        EQ2 = ( (2*dx)/(c *dT + dx) )  * ( Ez_Old(1,j) + Ez_Old(2,j) );
        EQ3 = ( (dx * (c * dT)^2 ) / (2 * dy * dy * (c * dT + dx)) ) * (MurCa + MurCb);
        Ez(1,j)   = -Ez_Old_Old(2,j) + EQ1 + EQ2 + EQ3;
    end
    
     for i = 2 : Nx - 1
        MurCa = Ez_Old(i+1,Ny) - 2 * Ez_Old(i,Ny) + Ez_Old(i-1,Ny);
        MurCb = Ez_Old(i+1,Ny-1) - 2 * Ez_Old(i,Ny-1) + Ez_Old(i-1,Ny-1);
        EQ1 = Mur * ( Ez(i,Ny-1) + Ez_Old_Old(i,Ny) );
        EQ2 = ( (2*dx)/(c *dT + dx) )  * ( Ez_Old(i,Ny) + Ez_Old(i,Ny-1) );
        EQ3 = ( (dx * (c * dT)^2 ) / (2 * dy * dy * (c * dT + dx)) ) * (MurCa + MurCb);
        Ez(i,Ny)   = -Ez_Old_Old(i,Ny-1) + EQ1 + EQ2 + EQ3;
     end
     
     for i = 2 : Nx - 1
        MurCa = Ez_Old(i+1,1) - 2 * Ez_Old(i,1) + Ez_Old(i-1,1);
        MurCb = Ez_Old(i+1,2) - 2 * Ez_Old(i,2) + Ez_Old(i-1,2);
        EQ1 = Mur * ( Ez(i,2) + Ez_Old_Old(i,1) );
        EQ2 = ( (2*dx)/(c *dT + dx) )  * ( Ez_Old(i,1) + Ez_Old(i,2) );
        EQ3 = ( (dx * (c * dT)^2 ) / (2 * dy * dy * (c * dT + dx)) ) * (MurCa + MurCb);
        Ez(i,1)   = -Ez_Old_Old(i,2) + EQ1 + EQ2 + EQ3;
     end
     
     for j = 2 : Ny - 1
        MurCa = Ez_Old(Nx,j+1) - 2 * Ez_Old(Nx,j) + Ez_Old(Nx,j-1);
        MurCb = Ez_Old(Nx-1,j+1) - 2 * Ez_Old(Nx-1,j) + Ez_Old(Nx-1,j-1);
        EQ1 = Mur * ( Ez(Nx-1,j) + Ez_Old_Old(Nx,j) );
        EQ2 = ( (2*dx)/(c *dT + dx) )  * ( Ez_Old(Nx,j) + Ez_Old(Nx-1,j) );
        EQ3 = ( (dx * (c * dT)^2 ) / (2 * dy * dy * (c * dT + dx)) ) * (MurCa + MurCb);
        Ez(Nx,j)   = -Ez_Old_Old(Nx-1,j) + EQ1 + EQ2 + EQ3;
     end

    for j = [1 Ny]
        Ez(1,j)   = Ez_Old(2,j) + Mur * (Ez(2,j) - Ez_Old(1,j));             
    end
    
    for i = [1 Nx]
        Ez(i,Ny)  = Ez_Old(i,Ny-1) + Mur * (Ez(i,Ny-1) - Ez_Old(i,Ny));             
    end
   
    for i = [1 Nx]
        Ez(i,1)   = Ez_Old(i,2) + Mur * (Ez(i,2) - Ez_Old(i,1));             
    end
    
    for j = [ 1 Ny]
       Ez(Nx,j)   = Ez_Old(Nx-1,j) + Mur * ( Ez(Nx-1,j) - Ez_Old(Nx, j));             
    end
  % ------------- END OF E-FIELD ABC -----------------------------------

     % Record Old Fields
   Ez_Old_Old_Old_Old = Ez_Old_Old_Old;
   Ez_Old_Old_Old = Ez_Old_Old;
   Ez_Old_Old = Ez_Old;
   Ez_Old = Ez;
   % -----------------
   Dz_Old_Old_Old = Dz_Old_Old;
   Dz_Old_Old = Dz_Old;
   Dz_Old = Dz;
end

ind = 1;
for i=1:Nx
	for j=1:Ny
		if (Slice(i,j) == 0)
	            Ez_cmpx(i,j) = 0;
	     	end
 	end
end
 
figure
M = abs(Ez_cmpx);
M = M./max(max(M));
pcolor(M) ; shading flat
t = sprintf('Field Intensity at %f GHz', Freq ./1e9);
title(t);
daspect([1 1 1])


