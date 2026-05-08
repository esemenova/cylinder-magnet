%% Magnetic field of a uniformly axially magnetized cylinder
%  Option A : analytical reduction of the surface integrals.
%
%  Geometry
%  --------
%  The cylinder is placed BELOW the XY plane.  Its top face lies at
%  z = 0 and its bottom face at z = -L (the cylinder occupies -L <= z <= 0).
%  Observation planes are at z = h above the top face.
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
%  The gradient tensor is obtained by central differences on the
%  analytical (quadrature-noiseless) field.
%
%  Validation: on-axis closed form (near & intermediate field)
%              point-dipole formula (far field).
% =========================================================================

clear; clc; close all;

%% --------------------- Parameters ---------------------------------------
mu0   = 4*pi*1e-7;
Br    = 1.0;
M0    = Br/mu0;
R     = 35e-3;
L     = 60e-3;

% --- Geometry: top face at z = 0, bottom face at z = -L
zTop  = 0;
zBot  = zTop - L;          % = -L
zCtr  = zTop - L/2;        % cylinder center, used by the dipole formula

% --- Observation distances ABOVE the top face
hList = [1e-3, 2e-3];      % 1 mm and 2 mm

% Gauss-Legendre nodes/weights for the radial integral on [0,R]
Nq          = 96;
[rGL, wGL]  = gaussLegendre1D(Nq, 0, R);
rGL = rGL(:);  wGL = wGL(:);

%% --------------------- 1. On-axis validation ----------------------------
zAx  = linspace(zTop + 1e-4, zTop + 20*L, 200).';
[~, ~, Bz_n] = Bcyl(zeros(size(zAx)), zeros(size(zAx)), zAx, R, L, M0, zTop, rGL, wGL);

% Exact on-axis closed form for cylinder occupying [zBot, zTop]
Bz_a = (mu0*M0/2) * ( (zAx - zBot)./sqrt(R^2 + (zAx - zBot).^2) - ...
                      (zAx - zTop)./sqrt(R^2 + (zAx - zTop).^2) );

% Dipole far-field (cylinder centred at zCtr)
m_dip = M0 * pi*R^2 * L;
Bz_d  = mu0*m_dip ./ (2*pi*abs(zAx - zCtr).^3);

figure('Name','On-axis validation','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
subplot(2,1,1);
loglog((zAx - zTop)/L, abs(Bz_n),'k-' ,'LineWidth',1.6); hold on;
loglog((zAx - zTop)/L, abs(Bz_a),'r--','LineWidth',1.2);
loglog((zAx - zTop)/L, abs(Bz_d),'b:' ,'LineWidth',1.4);
grid on; xlabel('(z - z_{top}) / L'); ylabel('|B_z|  [T]');
legend('analytical (1-D Gauss-Legendre)','on-axis closed form', ...
       'dipole far field','Location','southwest');
title('B_z above the cylinder axis');

subplot(2,1,2);
loglog((zAx - zTop)/L, abs(Bz_n - Bz_a)./abs(Bz_a),'r-','LineWidth',1.4); hold on;
loglog((zAx - zTop)/L, abs(Bz_n - Bz_d)./abs(Bz_d),'b-','LineWidth',1.4);
grid on; xlabel('(z - z_{top}) / L'); ylabel('relative error');
legend('vs on-axis analytic (exact)','vs dipole approximation', ...
       'Location','east');
title('Relative error of the numerical solution');

%% --------------------- 2. XY-plane field maps at z = h ------------------
xMax = R + 5e-3;
N    = 80;                       % grid resolution (NxN)
xv   = linspace(-xMax, xMax, N);
yv   = linspace(-xMax, xMax, N);
[X,Y] = meshgrid(xv, yv);

hd = 1e-6;                       % central-difference step
nH = numel(hList);
res(nH) = struct();              % preallocate result store

% --- compute fields and gradients for every requested height -------------
for ih = 1:nH
    h  = hList(ih);
    z0 = zTop + h;
    Z  = z0*ones(size(X));

    fprintf('\n=== h = %.1f mm ===\n', h*1e3);
    fprintf('B field on %dx%d grid ... ', N, N); tic;
    [Bx, By, Bz] = Bcyl(X, Y, Z, R, L, M0, zTop, rGL, wGL);
    fprintf('%.3f s\n', toc);

    fprintf('Gradient tensor ... '); tic;
    [Bx_xp, By_xp, Bz_xp] = Bcyl(X+hd, Y,    Z,    R, L, M0, zTop, rGL, wGL);
    [Bx_xm, By_xm, Bz_xm] = Bcyl(X-hd, Y,    Z,    R, L, M0, zTop, rGL, wGL);
    [Bx_yp, By_yp, Bz_yp] = Bcyl(X,    Y+hd, Z,    R, L, M0, zTop, rGL, wGL);
    [Bx_ym, By_ym, Bz_ym] = Bcyl(X,    Y-hd, Z,    R, L, M0, zTop, rGL, wGL);
    [~,     ~,     Bz_zp] = Bcyl(X,    Y,    Z+hd, R, L, M0, zTop, rGL, wGL);
    [~,     ~,     Bz_zm] = Bcyl(X,    Y,    Z-hd, R, L, M0, zTop, rGL, wGL);

    res(ih).h      = h;
    res(ih).Bx     = Bx;
    res(ih).By     = By;
    res(ih).Bz     = Bz;
    res(ih).dBx_dx = (Bx_xp - Bx_xm)/(2*hd);
    res(ih).dBx_dy = (Bx_yp - Bx_ym)/(2*hd);
    res(ih).dBy_dy = (By_yp - By_ym)/(2*hd);
    res(ih).dBz_dx = (Bz_xp - Bz_xm)/(2*hd);
    res(ih).dBz_dy = (Bz_yp - Bz_ym)/(2*hd);
    res(ih).dBz_dz = (Bz_zp - Bz_zm)/(2*hd);
    fprintf('%.3f s\n', toc);
end

Xc   = X*1e2;  Yc = Y*1e2;       % grid in cm for the maps
cLim = [-0.5 0.5];
th   = linspace(0,2*pi,200);

% ------------------------------------------------------------------ Fig 2
% B components for both heights side-by-side  (rows = component, cols = h)
figure('Name','B components','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
tiledlayout(3, nH, 'Padding','compact','TileSpacing','compact');

% Row 1 : Bx
for ih = 1:nH
    nexttile;
    pcolor(Xc, Yc, res(ih).Bx); shading interp;
    axis equal tight; colormap(jet); caxis(cLim);
    cb = colorbar; cb.Label.String = 'B_x (T)';
    xlabel('x (cm)'); ylabel('y (cm)');
    title(sprintf('h = %g mm', res(ih).h*1e3));
end
% Row 2 : By
for ih = 1:nH
    nexttile;
    pcolor(Xc, Yc, res(ih).By); shading interp;
    axis equal tight; colormap(jet); caxis(cLim);
    cb = colorbar; cb.Label.String = 'B_y (T)';
    xlabel('x (cm)'); ylabel('y (cm)');
end
% Row 3 : Bz as 3D surface
for ih = 1:nH
    nexttile;
    surf(Xc, Yc, res(ih).Bz, 'EdgeColor','none');
    colormap(jet); caxis(cLim);
    cb = colorbar; cb.Label.String = 'B_z (T)';
    xlabel('x (cm)'); ylabel('y (cm)'); zlabel('B_z (T)');
    view(-37.5, 30); axis tight;
end
sgtitle('B components','FontWeight','bold');

% ------------------------------------------------------------------ Fig 3
% In-plane field direction for both heights
figure('Name','In-plane field direction','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
tiledlayout(1, nH, 'Padding','compact','TileSpacing','compact');
s = 1:6:N;     % every 4th point
for ih = 1:nH
    nexttile;
    contourf(X*1e3, Y*1e3, res(ih).Bz, 20, 'LineStyle','none'); hold on;
    plot(R*cos(th)*1e3, R*sin(th)*1e3, 'w--','LineWidth',1.2);
    mg = hypot(res(ih).Bx, res(ih).By); mg(mg==0) = 1;
    
    quiver(X(s,s)*1e3, Y(s,s)*1e3, res(ih).Bx(s,s)./mg(s,s), res(ih).By(s,s)./mg(s,s), ...
       0.8, 'w','LineWidth',2);

    %quiver(X*1e3, Y*1e3, res(ih).Bx./mg, res(ih).By./mg, 1.8, 'k','LineWidth',1.0);
    axis equal tight; colorbar; colormap(jet);
    xlabel('x [mm]'); ylabel('y [mm]');
    title(sprintf('h = %g mm', res(ih).h*1e3));
end
sgtitle('B_z (color) and in-plane (B_x,B_y) direction','FontWeight','bold');

% ------------------------------------------------------------------ Fig 4
% Gradient components per height (5 components + div B), still one figure
% per distance because of panel count.
for ih = 1:nH
    figure('Name',sprintf('Gradient components, h = %g mm', res(ih).h*1e3),...
           'Color','w');
    set(gcf,'WindowStyle','docked','Color',[1 1 1]);
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    plotMap(X,Y,res(ih).dBx_dx,'\partial B_x/\partial x  [T/m]');
    plotMap(X,Y,res(ih).dBy_dy,'\partial B_y/\partial y  [T/m]');
    plotMap(X,Y,res(ih).dBz_dz,'\partial B_z/\partial z  [T/m]');
    plotMap(X,Y,res(ih).dBz_dx,'\partial B_z/\partial x  [T/m]');
    plotMap(X,Y,res(ih).dBz_dy,'\partial B_z/\partial y  [T/m]');
    plotMap(X,Y,res(ih).dBx_dx+res(ih).dBy_dy+res(ih).dBz_dz, ...
            'div B (should be ~ 0)  [T/m]');
    sgtitle(sprintf('Field gradient components, h = %g mm', res(ih).h*1e3));
end

% ------------------------------------------------------------------ Fig 5
% Gradient components of |B| for both heights (rows = h, cols = component)
figure('Name','Gradient of |B|','Color','w');
set(gcf,'WindowStyle','docked','Color',[1 1 1]);
tiledlayout(nH, 3, 'Padding','compact','TileSpacing','compact');
for ih = 1:nH
    Bxh = res(ih).Bx; Byh = res(ih).By; Bzh = res(ih).Bz;
    Bmag = sqrt(Bxh.^2 + Byh.^2 + Bzh.^2);
    dBmag_dx = (Bxh.*res(ih).dBx_dx + Byh.*res(ih).dBx_dy + Bzh.*res(ih).dBz_dx)./Bmag;
    dBmag_dy = (Bxh.*res(ih).dBx_dy + Byh.*res(ih).dBy_dy + Bzh.*res(ih).dBz_dy)./Bmag;
    dBmag_dz = (Bxh.*res(ih).dBz_dx + Byh.*res(ih).dBz_dy + Bzh.*res(ih).dBz_dz)./Bmag;

    htag = sprintf('h = %g mm', res(ih).h*1e3);
    plotMap(X,Y,dBmag_dx,['\partial |B| / \partial x  [T/m],  ' htag]);
    plotMap(X,Y,dBmag_dy,['\partial |B| / \partial y  [T/m],  ' htag]);
    plotMap(X,Y,dBmag_dz,['\partial |B| / \partial z  [T/m],  ' htag]);
end
sgtitle('Gradient of |B|','FontWeight','bold');

% =========================================================================
%                            Local functions
% =========================================================================
function [Bx, By, Bz] = Bcyl(x, y, z, R, L, M0, zTop, rGL, wGL)
% Closed-form (analytic) magnetic field of a uniformly axially magnetized
% cylinder.  Top face at z = zTop, bottom face at z = zTop - L.
% Vectorised over arbitrary x,y,z arrays.
mu0  = 4*pi*1e-7;
zBot = zTop - L;

% top and bottom face contributions (charge density +M0 and -M0)
[Brho_t, Bz_t] = faceField(x, y, z - zTop, R, +M0, rGL, wGL);
[Brho_b, Bz_b] = faceField(x, y, z - zBot, R, -M0, rGL, wGL);

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

% On-axis special case: Bx = By = 0 by symmetry; use closed-form B_z.
mask = rho < eps_r;
if any(mask(:))
    Bx(mask) = 0;
    By(mask) = 0;
    zh = z(mask);
    Bz(mask) = (mu0*M0/2) * ( (zh - zBot)./sqrt(R^2 + (zh - zBot).^2) - ...
                              (zh - zTop)./sqrt(R^2 + (zh - zTop).^2) );
end
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
