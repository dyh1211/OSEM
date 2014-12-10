function [B] = Inverse_Matrix_Builder(xc,yc,dx,dy,Probes,freq,Em,Efunc,SrcNo)

	Uo = 12.56e-7; % Freespace Permability
	Eo = 8.854e-12; % Freespace Permittivity
	cmpx = sqrt(-1); 
	M  = length(Probes(:,1)); % The  number of listerning probes
	N =  length(xc); % The number of segments in the solution domain
	k = sqrt( (2*pi*freq)^2 * Uo * Em*Eo - cmpx * 2 * pi * freq * Uo * 0); % The complex wave number
    Einc = Efunc;
    
	for m=1:M
    		x = Probes(m,1); 
    		y = Probes(m,2);
    		p = sqrt( (xc - x).^2 + (yc  - y).^2 );
    		a = sqrt( (dx*dy)/pi );
    	B(m,:) = (-cmpx * pi * k * 0.5) .* Einc .* a .* besselj(1,k*a) .* besselh(0,2,k*p');
 	end
end

  
