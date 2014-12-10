clc
close all
clear all

Eo = 8.854e-12;
Uo = 12.56e-7;
cmpx = sqrt(-1);

% -- These are the parameters for the simulation. 
Lx   = 1e-3 * 75;
Ly   = 1e-3 * 75;
dx   = 5e-3;
dy   = 5e-3;
dx2  = dx / 2;
dy2  = dy / 2;
Nx   = Lx / dx;
Ny   = Ly / dy;
Em   = 1 - 0 * cmpx;
c    = 3e8;
Freq = 5e9;

ind = 1;
for i = 1 : Nx
    for j = 1 : Ny
              Xc = dx2 +  (i-1) * dx - Lx / 2;
	      Yc = dy2 +  (j-1) * dy - Ly / 2;
	      R = sqrt(Xc * Xc + Yc * Yc);  
              xc(ind) = Xc;
              yc(ind) = Yc;
	      if (Xc > -0.02) && (Xc < 0.02) && (Yc > -0.02) && (Yc < 0.02)
  	            Er(ind,1) = 1.5 * Em - cmpx * 0.1;
	      else 
		   Er(ind, 1) = 1;
		end 
              ind = ind + 1;
    end
end

% Position of the listerning probes - a circular array
  Rx_Xradius = 0.05;
  Rx_Yradius = 0.05;
  Rx_angle   = 30;
  
  Rx_N  = 360 / Rx_angle; 
  angle = 0 : Rx_angle : 360 - Rx_angle;
  for n = 1 : Rx_N
      Probes(n,1) = Rx_Yradius * sin(angle(n) * pi / 180);     
      Probes(n,2) = Rx_Xradius * cos(angle(n) * pi / 180);
 
  end

%-Plot the Geometry----------------------------------- 
figure
subplot(2,2,1)
surf_from_scatter(xc,yc,real(Er)); 
title('Real(Er)') ; 
axis off
daspect([1 1 1])
subplot(2,2,2)
surf_from_scatter(xc,yc,imag(Er));
title('Imag(Er)')
daspect([1 1 1])
axis off
%----------------------------------------------------- 

% ----------------------------------------------------
fprintf('Running Forward Model (MoM) \n')
% Compute the Scattered Electric Fields at the Listerning Probes
% Each Set of Fields are stacked vertically for each frequency and antenna
 
for f=1:length(Freq)
       k = sqrt( (2*pi*Freq(f))^2 * Uo * Em*Eo  - cmpx * 2 * pi * Freq(f) * Uo * 0);
        for p=1:length(Probes(:,1))
            fprintf('Computing Fields From Antenna %i [%i]\n',p,length(Probes(:,1)));
            fflush(stdout);
            [Es(p).F,Et(p).F] = MoM(xc,yc,dx,dy,Er,k,Probes,p); 
        end
end

Ez=[]; 
for n=1:length(Es)
    Ez = [Ez ; Es(n).F];
end


% --------- <BEGIN> Inverse Procedure -------------------------------------

clear xc yc
% -------------- Readjust dx dy for the reconsutrction process ------
dx = 5e-3; 
dy = 5e-3;
dx2 = dx / 2;
dy2 = dy / 2;
Nx = Lx / dx;
Ny = Ly / dy;


%--------------------------------------------------------------------------
ind = 1;
for i=1:Nx
    for j=1:Ny
             % Determine the x,y coordinates of each cell centre
              xc(ind) = dx2 +  (i-1) * dx - (Lx/2); 
              yc(ind) = dy2 +  (j-1) * dy - (Ly/2);
              ind = ind + 1;
    end
end

xcR = xc;
ycR = yc;
clear Er_n
Er_n = ones(length(xcR),1);

N = length(xc);
Niter = 5;
for n = 1 : Niter
	fprintf('Computing Iteration %i [%i] \n',n,Niter);
	tic
    	clear Ez_n
    	ind = 1;
    	for f = 1 : length(Freq)
        	k = sqrt( (2*pi*Freq(f))^2 * Uo * Em*Eo  - cmpx * 2 * pi * Freq(f) * Uo * 0);
        	for p=1:length(Probes(:,1))
        	    fprintf('\t Calculating Fields at Antenna %i [%i]\n',p,length(Probes(:,1)));
        	    fflush(stdout);
        	    [Es,Ez_n(p).Et] = MoM(xc,yc,dx,dy,Er_n,k,Probes,p); 
        	end
    	end	

	B = [];
	% Construct 
	ind = 1;
	for f=1:length(Freq)
	    for p=1:length(Probes(:,1))
	        [Bn(p).F] = Inverse_Matrix_Builder(xc, yc, dx, dy, Probes, Freq(f),Em,Ez_n(p).Et,p);
	    end
	end
	
	for p=1:length(Probes)
	   B = [B ; Bn(p).F];
	end

	Lambda = 0; % Readjust depending on problem complexity
	L = length(xc);     % Determine the Length 
	N = length(B(1,:)); % Determine the number of unknowns 

	Breg = [B  ; Lambda * eye(N)];
	Ereg = [Ez ; zeros(N,1)];

	Treg =  pinv(Ereg) *  Breg;
	Er_n = (Treg + 1)'; % Calculate Reconstructed Permittivities with regularisation. 
end
% --------- <END> Inverse Procedure -------------------------------------

%-Plot the Image----------------------------------- 
subplot(2,2,3)
surf_from_scatter(xc,yc,real(Er_n)); 
title('Reconstructed Real(Er)') ; 
axis off
daspect([1 1 1])
subplot(2,2,4)
surf_from_scatter(xc,yc,imag(Er_n)); 
daspect([1 1 1])
axis off
title('Reconstructed Imag(Er)') ; 
%----------------------------------------------------- 



