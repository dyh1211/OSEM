function [Es,Et] =MoM(xc,yc,dx,dy,Er,k,Probes,SrcNo)
% xc - An array of cells centres (x coordinate)
% yc - An array of cells centres (y coordinate)
% Er - An array of permittivity values for each cell
% km - The wave number
% Probes - Position of listerning Probes
% Phi - Direction of Incoming incident wave
 
cmpx = sqrt(-1);

N = length(xc); % The size of the Matrix 
C = zeros(N,N); % Initilise the Coefficent Matrix

% Construct the [C] Matrix and the [b] Vector
for m=1:N
    x = xc(m);
    y = yc(m);
    Er_m = Er(m);  
    xdash = xc;
    ydash = yc;
    p = sqrt( (x - xdash).^2 + (y - ydash).^2 );
    a = (sqrt( (dx*dy)/pi));
    Er_n = Er;
    
    % Construct the [C] Matrix - using Matlab trickery. 
    C(m,:) = (cmpx * pi * k * a * 0.5) .*  besselj(1,k * a) .*(Er_n-1)  .* besselh(0,2,k .* p');
    C(m,m) = 1 + (Er_m - 1) * (cmpx/2) * (pi * k * a *  besselh(1,2,a*k) - 2 * cmpx);
end
pho = sqrt( (xc-Probes(SrcNo,1)).^2 + (yc - Probes(SrcNo,2)).^2 );
b(:,1) = cmpx * k * 0.25 * besselh(0,2,k .* pho');
    
% Size for the Total Field Inside the Domain
Et = C\b;  

% Initilise Vector to Store the Fields outside the Domain
Es = zeros(length(Probes(:,1)),1);

% Using the Fields inside the domain, compute the Fields at the listerning
% probes. By treating each cell as a point source we sum that up to find
% the contribution. 
 for n =1:length(Probes(:,1))
     x = Probes(n,1);
     y = Probes(n,2);
    
    total = 0;
    for i=1:length(b)
         Er_i = Er(i);
         xdash = xc(i);
         ydash = yc(i);
         p = sqrt( (x-xdash)^2 + (y - ydash)^2 );
         total= total + (-cmpx *pi*k/2) *  (Er_i -1) * Et(i) * a * besselj(1,k*a) * besselh(0,2,k * p);
     end
     % The scattered field at the n listerning probe. 
     Es(n,1) =  total ;
 end
 
