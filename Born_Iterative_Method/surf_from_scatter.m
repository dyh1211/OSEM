function surf_from_scatter(x,y,z)
	tri = delaunay(x,y);
	plot(x,y,'.')
	[r,c] = size(tri);
	h = trisurf(tri, x, y, z);
	axis off
	% Only Works in Matlab ------
	% axis vis3d
	% lighting phong;
	% ---------------------------
	shading interp;
	colorbar
	view(90,90);
end

