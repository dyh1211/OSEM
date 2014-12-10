close all
clear all
clc


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
	            Er(ind,1) = 40 * Em - cmpx * 0.1;
	      else 
		   Er(ind, 1) = 1;
	      end 
              ind = ind + 1;
    end
end

% Position of the Transmitter
Probe = [0 0.05];
  
%-Plot the Geometry----------------------------------- 
figure
subplot(2,2,1)
surf_from_scatter(xc, yc, real(Er)); 
title('Real(Er)') ; 
axis off
daspect([1 1 1])
subplot(2,2,2)
surf_from_scatter(xc, yc, imag(Er));
title('Imag(Er)')
daspect([1 1 1])
axis off
%----------------------------------------------------- 

% ----------------------------------------------------
 
for f=1:length(Freq)
       k = sqrt( (2*pi*Freq(f))^2 * Uo * Em*Eo  - cmpx * 2 * pi * Freq(f) * Uo * 0);
       fflush(stdout);
       [Es, Et] = MoM(xc, yc, dx, dy, Er, k, Probe, 1); 
end

subplot(2,2,3)
surf_from_scatter(xc,yc,real(Et));
title('Real(E) [V/m]')
daspect([1 1 1])
axis off

subplot(2,2,4)
surf_from_scatter(xc, yc, imag(Et));
title('Imag(E) [V/m]')
daspect([1 1 1])
axis off



