%% Magnetic field of a uniformly axially magnetized cylinder
%
%  Method:
%    Scalar magnetic potential phi_m as a surface integral over the two
%    flat faces of the cylinder (lateral surface contributes nothing for
%    axial magnetization).  The field H = -grad(phi_m) and the gradient
%    tensor dB_a/dx_b are obtained by differentiating analytically under
%    the integral sign and evaluating the resulting 2-D integrals with
%    MATLAB's adaptive integral2.
%
%  Validation:
%    - On-axis closed form  (near and intermediate field).
%    - Point-dipole formula (far field).
% =========================================================================

clear; clc; close all;

%% --------------------- Parameters ---------------------------------------
mu0 = 4*pi*1e-7;             % vacuum permeability [T m / A]
Br  = 1.0;                   % remanent flux density [T]
M0  = Br/mu0;                % magnetization magnitude [A/m]
R   = 35e-3;                  % cylinder radius   [m]
L   = 60e-3;                 % cylinder length   [m]
z0add = 1e-3;                % distance from the top [m]

%% --------------------- 1. On-axis validation ----------------------------
zAx  = linspace(L/2 + 1e-4, 20*L, 200).';
Bz_n = zeros(size(zAx));
for k = 1:numel(zAx)
    B        = Bfield(0, 0, zAx(k), R, L, M0);
    Bz_n(k)  = B(3);
end
Bz_a = Bz_axis_analytic(zAx, R, L, M0, mu0);
Bz_d = Bz_dipole_axis  (zAx, R, L, M0, mu0);

figure('Name','On-axis validation','Color','w');
set(gcf, 'WindowStyle', 'docked', 'Color', [1 1 1])

subplot(2,1,1);
loglog(zAx/L, abs(Bz_n),'k-' ,'LineWidth',1.6); hold on;
loglog(zAx/L, abs(Bz_a),'r--','LineWidth',1.2);
loglog(zAx/L, abs(Bz_d),'b:' ,'LineWidth',1.4);
grid on; xlabel('z / L'); ylabel('|B_z|  [T]');
legend('numerical (surface integral)','on-axis analytic','dipole far field', ...
       'Location','southwest');
title('B_z on the cylinder axis');

subplot(2,1,2);
loglog(zAx/L, abs(Bz_n - Bz_a)./abs(Bz_a),'r-','LineWidth',1.4); hold on;
loglog(zAx/L, abs(Bz_n - Bz_d)./abs(Bz_d),'b-','LineWidth',1.4);
grid on; xlabel('z / L'); ylabel('relative error');
legend('vs on-axis analytic (exact)','vs dipole approximation', ...
       'Location','east');
title('Relative error of the numerical solution');

%% --------------------- 2. XY-plane field map at z = z_0 -----------------
z0   = L/2 + z0add; %0.5*L;          % half a cylinder length above the top face
xMax = R+5e-3; % 3*R;
N    = 50;                   % grid resolution (NxN)

xv = linspace(-xMax, xMax, N);
yv = linspace(-xMax, xMax, N);
[X,Y] = meshgrid(xv, yv);

Bx = zeros(N); By = zeros(N); Bz = zeros(N);
dBx_dx = zeros(N); dBy_dy = zeros(N); dBz_dz = zeros(N);
dBz_dx = zeros(N); dBz_dy = zeros(N); dBx_dy = zeros(N);

fprintf('Computing %dx%d field grid ... ', N, N); tic;
for i = 1:N
    for j = 1:N
        Bv = Bfield(X(i,j), Y(i,j), z0, R, L, M0);
        Bx(i,j) = Bv(1); By(i,j) = Bv(2); Bz(i,j) = Bv(3);

        G = gradB(X(i,j), Y(i,j), z0, R, L, M0);
        dBx_dx(i,j) = G(1,1);
        dBy_dy(i,j) = G(2,2);
        dBz_dz(i,j) = G(3,3);
        dBz_dx(i,j) = G(3,1);
        dBz_dy(i,j) = G(3,2);
        dBx_dy(i,j) = G(1,2);
    end
end
fprintf('%.1f s\n', toc);

%% --------------------- 3. Plots -----------------------------------------
% Figure A : B components -- Bx, By as 2D maps, Bz as 3D surface
figure('Name','B components','Color','w');
set(gcf, 'WindowStyle', 'docked', 'Color', [1 1 1])

Xc = X*1e2;  Yc = Y*1e2;            % grid in cm
cLim = [-0.5 0.5];                  % common color range

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

sgtitle(sprintf('h = %g mm', z0add*1e3), 'FontWeight','bold');

% Figure B : in-plane direction (quiver) over B_z color map
figure('Name','In-plane field direction','Color','w');
set(gcf, 'WindowStyle', 'docked', 'Color', [1 1 1])

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
set(gcf, 'WindowStyle', 'docked', 'Color', [1 1 1])

tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
plotMap(X,Y,dBx_dx,'\partial B_x/\partial x  [T/m]');
plotMap(X,Y,dBy_dy,'\partial B_y/\partial y  [T/m]');
plotMap(X,Y,dBz_dz,'\partial B_z/\partial z  [T/m]');
plotMap(X,Y,dBz_dx,'\partial B_z/\partial x  [T/m]');
plotMap(X,Y,dBz_dy,'\partial B_z/\partial y  [T/m]');
plotMap(X,Y,dBx_dx+dBy_dy+dBz_dz,'div B (should be ~ 0)  [T/m]');
sgtitle(sprintf('Field gradient components at z = %.2f mm', z0*1e3));

% Figure D : gradient components of |B|
% |B| = sqrt(Bx^2+By^2+Bz^2);  d|B|/da = (Bx*dBx/da + By*dBy/da + Bz*dBz/da)/|B|
% Using symmetry of the gradient tensor in source-free region:
%   dBy/dx = dBx/dy,  dBx/dz = dBz/dx,  dBy/dz = dBz/dy
Bmag    = sqrt(Bx.^2 + By.^2 + Bz.^2);
dBmag_dx = (Bx.*dBx_dx + By.*dBx_dy + Bz.*dBz_dx) ./ Bmag;
dBmag_dy = (Bx.*dBx_dy + By.*dBy_dy + Bz.*dBz_dy) ./ Bmag;
dBmag_dz = (Bx.*dBz_dx + By.*dBz_dy + Bz.*dBz_dz) ./ Bmag;

figure('Name','Gradient of |B|','Color','w');
set(gcf, 'WindowStyle', 'docked', 'Color', [1 1 1])
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
plotMap(X,Y,dBmag_dx,'\partial |B| / \partial x  [T/m]');
plotMap(X,Y,dBmag_dy,'\partial |B| / \partial y  [T/m]');
plotMap(X,Y,dBmag_dz,'\partial |B| / \partial z  [T/m]');
sgtitle(sprintf('Gradient of |B| at z = %.2f mm', z0*1e3));

% =========================================================================
%                            Local functions
% =========================================================================
function B = Bfield(x, y, z, R, L, M0)
% B = mu0 * H at (x,y,z) outside the magnet, from surface integral.
mu0 = 4*pi*1e-7;
[Hxt, Hyt, Hzt] = faceH(x, y, z, +L/2, +M0, R);
[Hxb, Hyb, Hzb] = faceH(x, y, z, -L/2, -M0, R);
B = mu0 * [Hxt+Hxb; Hyt+Hyb; Hzt+Hzb];
end

function [Hx, Hy, Hz] = faceH(x, y, z, zp, sigma, R)
% Field contribution of a uniformly charged disc at z = zp.
%   H_a = (sigma/4pi) * iint (a - a')/d^3 * r' dr' dphi'
fX = @(rp,php) (x - rp.*cos(php)) ./ ...
       ((x - rp.*cos(php)).^2 + (y - rp.*sin(php)).^2 + (z - zp).^2).^1.5 .* rp;
fY = @(rp,php) (y - rp.*sin(php)) ./ ...
       ((x - rp.*cos(php)).^2 + (y - rp.*sin(php)).^2 + (z - zp).^2).^1.5 .* rp;
fZ = @(rp,php) (z - zp) ./ ...
       ((x - rp.*cos(php)).^2 + (y - rp.*sin(php)).^2 + (z - zp).^2).^1.5 .* rp;
opts = {'AbsTol',1e-10,'RelTol',1e-7};
c    = sigma/(4*pi);
Hx = c * integral2(fX, 0, R, 0, 2*pi, opts{:});
Hy = c * integral2(fY, 0, R, 0, 2*pi, opts{:});
Hz = c * integral2(fZ, 0, R, 0, 2*pi, opts{:});
end

function G = gradB(x, y, z, R, L, M0)
% 3x3 gradient tensor  G(a,b) = dB_a / dx_b  at (x,y,z).
mu0 = 4*pi*1e-7;
Gt  = faceGradH(x, y, z, +L/2, +M0, R);
Gb  = faceGradH(x, y, z, -L/2, -M0, R);
G   = mu0 * (Gt + Gb);
end

function G = faceGradH(x, y, z, zp, sigma, R)
% dH_a/dx_b = (sigma/4pi) * iint [ -delta_ab/d^3 + 3*(a-a')(b-b')/d^5 ] r' dr' dphi'
ax = @(rp,php) x - rp.*cos(php);
ay = @(rp,php) y - rp.*sin(php);
az = @(rp,php) (z - zp) + 0.*rp;
d  = @(rp,php) sqrt(ax(rp,php).^2 + ay(rp,php).^2 + az(rp,php).^2);
co = {ax, ay, az};

opts = {'AbsTol',1e-8,'RelTol',1e-6};
c    = sigma/(4*pi);
G    = zeros(3,3);
for a = 1:3
    for b = a:3
        if a == b
            ker = @(rp,php) ( -1./d(rp,php).^3 + ...
                              3*co{a}(rp,php).^2 ./ d(rp,php).^5 ) .* rp;
        else
            ker = @(rp,php) ( 3*co{a}(rp,php).*co{b}(rp,php) ./ ...
                              d(rp,php).^5 ) .* rp;
        end
        v = c * integral2(ker, 0, R, 0, 2*pi, opts{:});
        G(a,b) = v;
        G(b,a) = v;     % gradient tensor is symmetric outside sources
    end
end
end

function Bz = Bz_axis_analytic(z, R, L, M0, mu0)
% Exact on-axis closed form (rho = 0).
zp = z + L/2;  zm = z - L/2;
Bz = (mu0*M0/2) * ( zp ./ sqrt(R^2 + zp.^2) - zm ./ sqrt(R^2 + zm.^2) );
end

function Bz = Bz_dipole_axis(z, R, L, M0, mu0)
% Far-field dipole approximation along the axis.
m  = M0 * pi*R^2 * L;                 % total dipole moment
Bz = mu0 * m ./ (2*pi * abs(z).^3);   % B_z on axis of a z-dipole
end

function plotMap(X, Y, F, ttl)
nexttile;
contourf(X*1e3, Y*1e3, F, 20, 'LineStyle','none');
colormap(jet)
axis equal tight; colorbar;
xlabel('x [mm]'); ylabel('y [mm]'); title(ttl);
end
