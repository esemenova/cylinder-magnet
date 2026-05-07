%% Magnetic field of a uniformly axially magnetized cylinder
%  Option A : analytical reduction of the surface integrals.
%
%  Method
%  ------
%  Each end face of the cylinder is a uniformly charged disc.  The
%  magnetic-field components H_alpha = (sigma/4pi) * iint (r-r')/d^3 dS'
%  are derived analytically *under the integral sign*:
%
%    1) The angular integral over phi' is performed in closed form
%       and yields complete elliptic integrals K(m) and E(m) with
%       m(r',rho,zeta) = 4*r'*rho / ((r'+rho)^2 + zeta^2).
%    2) The remaining 1-D integral over r' on [0,R] is evaluated by
%       Gauss-Legendre quadrature, fully vectorised over the whole
%       observation grid.
%
%  This replaces the per-point integral2() of the brute-force script
%  by a single matrix product, giving a speed-up of ~10^3.
%
%  Field formulas used (derived from the surface integral):
%
%    H_z(rho,zeta) = (sigma/pi) * zeta *
%                    int_0^R r' E(m) / [ ((r'-rho)^2+zeta^2) *
%                                        sqrt((r'+rho)^2+zeta^2) ] dr'
%
%    H_rho(rho,zeta) = sigma/(2*pi*rho) *
%        int_0^R r' / [ ((r'-rho)^2+zeta^2)*sqrt((r'+rho)^2+zeta^2) ] *
%             [ (rho^2 - r'^2 - zeta^2) E(m)
%             + ((r'-rho)^2 + zeta^2)  K(m) ] dr'
%
%  with m = 4 r' rho / ((r'+rho)^2 + zeta^2).
%
%  The cylinder field is the sum of the two face contributions
%   (top  face : z' = +L/2,  sigma = +M0
%    bot  face : z' = -L/2,  sigma = -M0).
%
%  The gradient tensor d B_a / d x_b is obtained by central differences
%  taken on the analytical (quadrature-noiseless) field with a small step
%  h ~ 1e-6 m, accuracy ~1e-9.
%
%  Validation:
%   - on-axis closed form  (near & intermediate field)
%   - point-dipole formula (far field)
% =========================================================================

clear; clc; close all;

%% --------------------- Parameters ---------------------------------------
mu0   = 4*pi*1e-7;
Br    = 1.0;
M0    = Br/mu0;
R     = 35e-3;
L     = 60e-3;
b     = L/2;
z0add = 1e-3;

% Gauss-Legendre nodes/weights for the radial integral on [0,R]
Nq          = 96;
[rGL, wGL]  = gaussLegendre1D(Nq, 0, R);
rGL = rGL(:);  wGL = wGL(:);

%% --------------------- 1. On-axis validation ----------------------------
zAx  = linspace(b+1e-4, 20*L, 200).';
[~, ~, Bz_n] = Bcyl(zeros(size(zAx)), zeros(size(zAx)), zAx, R, b, M0, rGL, wGL);

Bz_a = (mu0*M0/2) * ( (zAx+b)./sqrt(R^2+(zAx+b).^2) - ...
                      (zAx-b)./sqrt(R^2+(zAx-b).^2) );
m_dip = M0 * pi*R^2 * L;
Bz_d  = mu0*m_dip ./ (2*pi*abs(zAx).^3);

figure('Name','On-axis validation','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
subplot(2,1,1);
loglog(zAx/L, abs(Bz_n),'k-' ,'LineWidth',1.6); hold on;
loglog(zAx/L, abs(Bz_a),'r--','LineWidth',1.2);
loglog(zAx/L, abs(Bz_d),'b:' ,'LineWidth',1.4);
grid on; xlabel('z / L'); ylabel('|B_z|  [T]');
legend('analytical (1-D Gauss-Legendre)','on-axis closed form', ...
       'dipole far field','Location','southwest');
title('B_z on the cylinder axis');

subplot(2,1,2);
loglog(zAx/L, abs(Bz_n - Bz_a)./abs(Bz_a),'r-','LineWidth',1.4); hold on;
loglog(zAx/L, abs(Bz_n - Bz_d)./abs(Bz_d),'b-','LineWidth',1.4);
grid on; xlabel('z / L'); ylabel('relative error');
legend('vs on-axis analytic (exact)','vs dipole approximation', ...
       'Location','east');
title('Relative error of the numerical solution');

%% --------------------- 2. XY-plane field map at z = z_0 -----------------
z0   = b + z0add;
xMax = R + 5e-3;
N    = 80;                       % grid resolution (NxN)

xv = linspace(-xMax, xMax, N);
yv = linspace(-xMax, xMax, N);
[X,Y] = meshgrid(xv, yv);
Z     = z0*ones(size(X));

fprintf('B field on %dx%d grid ... ', N, N); tic;
[Bx, By, Bz] = Bcyl(X, Y, Z, R, b, M0, rGL, wGL);
fprintf('%.3f s\n', toc);

%% Gradient tensor : central differences on the (analytic) B
hd = 1e-6;
fprintf('Gradient tensor ... '); tic;
[Bx_xp, By_xp, Bz_xp] = Bcyl(X+hd, Y,    Z,    R, b, M0, rGL, wGL);
[Bx_xm, By_xm, Bz_xm] = Bcyl(X-hd, Y,    Z,    R, b, M0, rGL, wGL);
[Bx_yp, By_yp, Bz_yp] = Bcyl(X,    Y+hd, Z,    R, b, M0, rGL, wGL);
[Bx_ym, By_ym, Bz_ym] = Bcyl(X,    Y-hd, Z,    R, b, M0, rGL, wGL);
[~,     ~,     Bz_zp] = Bcyl(X,    Y,    Z+hd, R, b, M0, rGL, wGL);
[~,     ~,     Bz_zm] = Bcyl(X,    Y,    Z-hd, R, b, M0, rGL, wGL);

dBx_dx = (Bx_xp - Bx_xm)/(2*hd);
dBx_dy = (Bx_yp - Bx_ym)/(2*hd);
dBy_dy = (By_yp - By_ym)/(2*hd);
dBz_dx = (Bz_xp - Bz_xm)/(2*hd);
dBz_dy = (Bz_yp - Bz_ym)/(2*hd);
dBz_dz = (Bz_zp - Bz_zm)/(2*hd);
fprintf('%.3f s\n', toc);

%% --------------------- 3. Plots -----------------------------------------
% Figure A : Bx, By as 2D maps, Bz as 3D surface
figure('Name','B components','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);

Xc = X*1e2; Yc = Y*1e2;
cLim = [-0.5 0.5];

tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

nexttile;
pcolor(Xc, Yc, Bx); shading interp;
axis equal tight; colormap(jet); caxis(cLim);
cb = colorbar; cb.Label.String = 'B_x (T)';
xlabel('x (cm)'); ylabel('y (cm)');

nexttile;
pcolor(Xc, Yc, By); shading interp;
axis equal tight; colormap(jet); caxis(cLim);
cb = colorbar; cb.Label.String = 'B_y (T)';
xlabel('x (cm)'); ylabel('y (cm)');

nexttile;
surf(Xc, Yc, Bz, 'EdgeColor','none');
colormap(jet); caxis(cLim);
cb = colorbar; cb.Label.String = 'B_z (T)';
xlabel('x (cm)'); ylabel('y (cm)'); zlabel('B_z (T)');
view(-37.5, 30); axis tight;

sgtitle(sprintf('h = %g mm', z0add*1e3),'FontWeight','bold');

% Figure B : in-plane field direction
figure('Name','In-plane field direction','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
contourf(X*1e3, Y*1e3, Bz, 20, 'LineStyle','none'); hold on;
th = linspace(0,2*pi,200);
plot(R*cos(th)*1e3, R*sin(th)*1e3, 'w--','LineWidth',1.2);
mag = hypot(Bx,By); mag(mag==0) = 1;
quiver(X*1e3, Y*1e3, Bx./mag, By./mag, 0.5, 'k', 'LineWidth', 0.7);
axis equal tight; colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('B_z (color) and (B_x,B_y) direction, z = %.2f mm', z0*1e3));

% Figure C : gradient components and divergence check
figure('Name','Gradient components','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
plotMap(X,Y,dBx_dx,'\partial B_x/\partial x  [T/m]');
plotMap(X,Y,dBy_dy,'\partial B_y/\partial y  [T/m]');
plotMap(X,Y,dBz_dz,'\partial B_z/\partial z  [T/m]');
plotMap(X,Y,dBz_dx,'\partial B_z/\partial x  [T/m]');
plotMap(X,Y,dBz_dy,'\partial B_z/\partial y  [T/m]');
plotMap(X,Y,dBx_dx+dBy_dy+dBz_dz,'div B (should be ~ 0)  [T/m]');
sgtitle(sprintf('Field gradient components at z = %.2f mm', z0*1e3));

% Figure D : gradient components of |B|
Bmag    = sqrt(Bx.^2 + By.^2 + Bz.^2);
dBmag_dx = (Bx.*dBx_dx + By.*dBx_dy + Bz.*dBz_dx) ./ Bmag;
dBmag_dy = (Bx.*dBx_dy + By.*dBy_dy + Bz.*dBz_dy) ./ Bmag;
dBmag_dz = (Bx.*dBz_dx + By.*dBz_dy + Bz.*dBz_dz) ./ Bmag;
figure('Name','Gradient of |B|','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
plotMap(X,Y,dBmag_dx,'\partial |B| / \partial x  [T/m]');
plotMap(X,Y,dBmag_dy,'\partial |B| / \partial y  [T/m]');
plotMap(X,Y,dBmag_dz,'\partial |B| / \partial z  [T/m]');
sgtitle(sprintf('Gradient of |B| at z = %.2f mm', z0*1e3));

% =========================================================================
%                            Local functions
% =========================================================================
function [Bx, By, Bz] = Bcyl(x, y, z, R, b, M0, rGL, wGL)
% Closed-form (analytic) magnetic field of a uniformly axially magnetized
% cylinder, vectorised over arbitrary x,y,z arrays.
mu0 = 4*pi*1e-7;

% top and bottom face contributions (charge density M0 and -M0)
[Brho_t, Bz_t] = faceField(x, y, z - b, R, +M0, rGL, wGL);
[Brho_b, Bz_b] = faceField(x, y, z + b, R, -M0, rGL, wGL);

Brho = mu0 * (Brho_t + Brho_b);
Bz   = mu0 * (Bz_t   + Bz_b  );

rho   = sqrt(x.^2 + y.^2);
eps_r = 1e-12 * R;

% Cartesian projection (Bphi = 0 by axial symmetry)
cosphi = zeros(size(rho));
sinphi = zeros(size(rho));
nz     = rho > eps_r;
cosphi(nz) = x(nz) ./ rho(nz);
sinphi(nz) = y(nz) ./ rho(nz);
Bx = Brho .* cosphi;
By = Brho .* sinphi;
end

function [Brho, Bz] = faceField(x, y, zeta, R, sigma, rGL, wGL)
% Field of one charged disc of radius R, surface charge sigma, at axial
% offset zeta from the observation plane.  Phi-integration done analytically
% (gives K,E); radial integration by Gauss-Legendre quadrature.
sz = size(x);
rho   = sqrt(x.^2 + y.^2);
rho_v = rho(:).';                 % 1 x Nobs
zet_v = zeta(:).';                % 1 x Nobs

% Auxiliary arrays  Nq x Nobs
Aplus  = (rGL + rho_v).^2 + zet_v.^2;       % (r'+rho)^2 + zeta^2
Aminus = (rGL - rho_v).^2 + zet_v.^2;       % (r'-rho)^2 + zeta^2
m      = 4*rGL.*rho_v ./ Aplus;             % elliptic parameter k^2
m      = min(max(m, 0), 1 - 1e-14);
[K, E] = ellipke(m);

denom  = Aminus .* sqrt(Aplus);             % nonzero outside the magnet

% Integrand for H_z
intHz   = rGL .* E ./ denom;

% Integrand for H_rho  ( the 1/rho is factored out in front )
bracket = (rho_v.^2 - rGL.^2 - zet_v.^2).*E + Aminus.*K;
intHrho = (rGL ./ denom) .* bracket;

% Quadrature (sum over Nq with weights wGL)
Hz_v   = (sigma/pi)        * zet_v                 .* sum(wGL .* intHz  , 1);
Hrho_v = (sigma/(2*pi))   ./ max(rho_v, 1e-30)    .* sum(wGL .* intHrho, 1);

% Force H_rho = 0 on the axis
Hrho_v(rho_v < 1e-12*R) = 0;

Brho = reshape(Hrho_v, sz);
Bz   = reshape(Hz_v  , sz);
end

function [x, w] = gaussLegendre1D(N, a, b)
% Golub-Welsch: N-point Gauss-Legendre nodes and weights on [a,b].
beta = .5 ./ sqrt(1 - (2*(1:N-1)).^(-2));
T    = diag(beta,1) + diag(beta,-1);
[V, D] = eig(T);
[x, idx] = sort(diag(D));
w  = 2 * V(1, idx).^2;
x  = 0.5*(b-a)*x + 0.5*(b+a);
w  = 0.5*(b-a)*w;
end

function plotMap(X, Y, F, ttl)
nexttile;
contourf(X*1e3, Y*1e3, F, 20, 'LineStyle','none');
axis equal tight; colorbar; colormap(jet)
xlabel('x [mm]'); ylabel('y [mm]'); title(ttl);
end
