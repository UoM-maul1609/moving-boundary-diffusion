% Reproduce the moving-boundary diffusion example from output.nc.
% The default namelist.in is configured for this demonstration.

ncfile = '/tmp/output.nc';

time   = double(ncread(ncfile,'time'));
r      = double(ncread(ncfile,'r'));
c      = double(ncread(ncfile,'c'));
radius = double(ncread(ncfile,'radius'));

% Files written by diffusion.f90 are [kp,ncomp,time] for c and [kp,time] for r.
water = squeeze(c(:,1,:));
solute = squeeze(c(:,2,:));
total = water + solute;
xwater = water ./ total;
xwater(total <= 0) = NaN;

figure;
pcolor(time, r, xwater);
shading flat;
colormap(parula);
cb = colorbar;
ylabel(cb,'mole fraction');
xlabel('seconds (s)');
ylabel('radius (m)');
ylim([0 3.e-7]);
box on;

% Optional consistency check: the scalar boundary should equal the active
% moving face throughout growth and evaporation after the shrink-radius fix.
fprintf('radius range: %.6g to %.6g m\n', min(radius), max(radius));
