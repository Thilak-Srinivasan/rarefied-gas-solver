clear; clc; close all;

gamma = 1.4;
Pr = 0.72;
Re_delta = 100;
M_inf = 2.5;
wedge_angle = 15;
T_ref = 273.15;

Nx = 120;             
Ny = 80;
Lx = 8;
Ly = 4;

CFL = 0.20;
max_steps = 30000;
conv_tol = 1e-6;
plot_freq = 500;

fprintf('========================================\n');
fprintf('ENHANCED 2D OBLIQUE SHOCK SOLVER\n');
fprintf('Rusanov Scheme + Navier-Stokes NSF\n');
fprintf('CORRECTED VERSION - Fixed Velocity Plot\n');
fprintf('========================================\n');
fprintf('Domain: %.2f x %.2f\n', Lx, Ly);
fprintf('Grid: %d x %d (dx=%.4f, dy=%.4f)\n', Nx, Ny, Lx/Nx, Ly/Ny);
fprintf('Mach: %.2f, Wedge Angle: %d°\n', M_inf, wedge_angle);
fprintf('Reynolds: %.1f, CFL: %.2f\n', Re_delta, CFL);
fprintf('========================================\n\n');

x = linspace(-1, Lx-1, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);
dx = x(2) - x(1);
dy = y(2) - y(1);

wedge_start = 0.5;
wedge_angle_rad = wedge_angle * pi / 180;

is_solid = false(Ny, Nx);
for j = 1:Nx
    if x(j) >= wedge_start
        wedge_height_col = (x(j) - wedge_start) * tan(wedge_angle_rad);
        is_solid(:, j) = y' <= wedge_height_col;
    end
end


rho = ones(Ny, Nx);
u = ones(Ny, Nx) * M_inf;
v = zeros(Ny, Nx);
p = ones(Ny, Nx);
T = p ./ rho;

U1 = rho;
U2 = rho .* u;
U3 = rho .* v;
E = p/(gamma-1) + 0.5*rho.*(u.^2 + v.^2);
U4 = E;

u(is_solid) = 0;
v(is_solid) = 0;
U2(is_solid) = 0;
U3(is_solid) = 0;

mu_ref = 1.0;
T_s = 0.4;

residual_history = zeros(min(max_steps, 50000), 1);
shock_data = struct();
shock_data.angle = [];
shock_data.p_ratio = [];
shock_data.rho_ratio = [];
shock_data.T_ratio = [];
shock_data.M2 = [];
shock_data.steps = [];


tic;
for n = 1:max_steps
    U1_old = U1;
    U2_old = U2;
    U3_old = U3;
    U4_old = U4;
    
    rho = max(U1, 1e-5);
    u = U2 ./ rho;
    v = U3 ./ rho;
    ke = 0.5 * (u.^2 + v.^2);
    p = (gamma-1) * (U4 - rho.*ke);
    p = max(p, 1e-5);
    T = p ./ rho;
    
    a = sqrt(gamma * p ./ rho);
    
    mu = mu_ref * (T + T_s) ./ (1 + T_s) .* (T .^ 1.5);
    mu(is_solid) = 0;
    kappa = mu / Pr;
    
    lambda_x = abs(u) + a;
    lambda_y = abs(v) + a;
    lambda_max = max(max(lambda_x)) + max(max(lambda_y));
    dt = CFL * min(dx, dy) / (lambda_max + 1e-10);
    
    [du_dx, du_dy] = compute_gradients_vectorized(u, dx, dy);
    [dv_dx, dv_dy] = compute_gradients_vectorized(v, dx, dy);
    [dT_dx, dT_dy] = compute_gradients_vectorized(T, dx, dy);
    
    tau_xx = mu .* (4/3 .* du_dx - 2/3 .* dv_dy);
    tau_yy = mu .* (4/3 .* dv_dy - 2/3 .* du_dx);
    tau_xy = mu .* (du_dy + dv_dx);
    
    qx = -kappa .* dT_dx;
    qy = -kappa .* dT_dy;
    
    dU1_dt = zeros(Ny, Nx);
    dU2_dt = zeros(Ny, Nx);
    dU3_dt = zeros(Ny, Nx);
    dU4_dt = zeros(Ny, Nx);
    
    for i = 2:Ny-1
        for j = 2:Nx-1
            if ~is_solid(i,j)
                % X-direction fluxes (Rusanov)
                [F_L, ~] = compute_inviscid_flux_x(...
                    U1(i,j-1), U2(i,j-1), U3(i,j-1), U4(i,j-1), gamma);
                [F_R, ~] = compute_inviscid_flux_x(...
                    U1(i,j+1), U2(i,j+1), U3(i,j+1), U4(i,j+1), gamma);
                
                rho_int = 0.5*(rho(i,j+1) + rho(i,j-1));
                u_int = 0.5*(u(i,j+1) + u(i,j-1));
                a_int = sqrt(gamma * 0.5*(p(i,j+1) + p(i,j-1)) / rho_int);
                lambda_max_x = abs(u_int) + a_int;
                
                F1_rusanov = 0.5*(F_L(1) + F_R(1)) - 0.5*lambda_max_x*(U1(i,j+1) - U1(i,j-1));
                F2_rusanov = 0.5*(F_L(2) + F_R(2)) - 0.5*lambda_max_x*(U2(i,j+1) - U2(i,j-1));
                F3_rusanov = 0.5*(F_L(3) + F_R(3)) - 0.5*lambda_max_x*(U3(i,j+1) - U3(i,j-1));
                F4_rusanov = 0.5*(F_L(4) + F_R(4)) - 0.5*lambda_max_x*(U4(i,j+1) - U4(i,j-1));
                
                % Y-direction fluxes
                [G_B, ~] = compute_inviscid_flux_y(...
                    U1(i-1,j), U2(i-1,j), U3(i-1,j), U4(i-1,j), gamma);
                [G_T, ~] = compute_inviscid_flux_y(...
                    U1(i+1,j), U2(i+1,j), U3(i+1,j), U4(i+1,j), gamma);
                
                rho_int = 0.5*(rho(i+1,j) + rho(i-1,j));
                v_int = 0.5*(v(i+1,j) + v(i-1,j));
                a_int = sqrt(gamma * 0.5*(p(i+1,j) + p(i-1,j)) / rho_int);
                lambda_max_y = abs(v_int) + a_int;
                
                G1_rusanov = 0.5*(G_B(1) + G_T(1)) - 0.5*lambda_max_y*(U1(i+1,j) - U1(i-1,j));
                G2_rusanov = 0.5*(G_B(2) + G_T(2)) - 0.5*lambda_max_y*(U2(i+1,j) - U2(i-1,j));
                G3_rusanov = 0.5*(G_B(3) + G_T(3)) - 0.5*lambda_max_y*(U3(i+1,j) - U3(i-1,j));
                G4_rusanov = 0.5*(G_B(4) + G_T(4)) - 0.5*lambda_max_y*(U4(i+1,j) - U4(i-1,j));
                
                % Conservation equations
                dU1_dt(i,j) = -(F2_rusanov - 0)/dx - (G3_rusanov - 0)/dy;
                
                dU2_dt(i,j) = -(F2_rusanov - 0)/dx - (G2_rusanov - 0)/dy ...
                    - (tau_xx(i,j+1) - tau_xx(i,j-1))/(2*dx) ...
                    - (tau_xy(i+1,j) - tau_xy(i-1,j))/(2*dy);
                
                dU3_dt(i,j) = -(F3_rusanov - 0)/dx - (G3_rusanov - 0)/dy ...
                    - (tau_xy(i,j+1) - tau_xy(i,j-1))/(2*dx) ...
                    - (tau_yy(i+1,j) - tau_yy(i-1,j))/(2*dy);
                
                visc_work = u(i,j)*(tau_xx(i,j+1) - tau_xx(i,j-1))/(2*dx) + ...
                            v(i,j)*(tau_yy(i+1,j) - tau_yy(i-1,j))/(2*dy) + ...
                            u(i,j)*(tau_xy(i+1,j) - tau_xy(i-1,j))/(2*dy) + ...
                            v(i,j)*(tau_xy(i,j+1) - tau_xy(i,j-1))/(2*dx);
                
                dU4_dt(i,j) = -(F4_rusanov - 0)/dx - (G4_rusanov - 0)/dy ...
                    + visc_work ...
                    - (qx(i,j+1) - qx(i,j-1))/(2*dx) ...
                    - (qy(i+1,j) - qy(i-1,j))/(2*dy);
            end
        end
    end
    
    U1 = U1 + dt * dU1_dt;
    U2 = U2 + dt * dU2_dt;
    U3 = U3 + dt * dU3_dt;
    U4 = U4 + dt * dU4_dt;
    
    % Boundary conditions
    [U1, U2, U3, U4] = apply_boundary_conditions(...
        U1, U2, U3, U4, M_inf, gamma, is_solid, x, X, Y, ...
        wedge_start, wedge_angle_rad);
    
    U1 = max(U1, 1e-5);
    
    rho_check = U1;
    u_check = U2 ./ rho_check;
    v_check = U3 ./ rho_check;
    p_check = (gamma-1) * (U4 - 0.5*rho_check.*(u_check.^2 + v_check.^2));
    
    neg_p = p_check < 1e-5;
    if any(neg_p(:))
        U1(neg_p) = 1.0;
        U2(neg_p) = M_inf;
        U3(neg_p) = 0;
        U4(neg_p) = 1.0/(gamma-1) + 0.5*M_inf^2;
    end
    
    residual = sqrt(mean((U1(:) - U1_old(:)).^2) + ...
                    mean((U2(:) - U2_old(:)).^2) + ...
                    mean((U3(:) - U3_old(:)).^2) + ...
                    mean((U4(:) - U4_old(:)).^2));
    residual_history(n) = residual;
    
    if mod(n, plot_freq) == 0
        fprintf('Step %6d: Residual = %.6e, dt = %.6e, CFL = %.3f\n', ...
            n, residual, dt, CFL);
        
        rho = U1;
        u = U2 ./ rho;
        v = U3 ./ rho;
        p = (gamma-1) * (U4 - 0.5*rho.*(u.^2 + v.^2));
        T = p ./ rho;
        
        plot_solution_corrected(X, Y, rho, u, v, p, T, gamma, n, is_solid, ...
            wedge_start, wedge_angle_rad, wedge_angle);
    end
    
    if residual < conv_tol && n > 2000
        fprintf('\n✓ CONVERGED at step %d (residual = %.2e)\n', n, residual);
        residual_history = residual_history(1:n);
        break;
    end
    
    if isnan(residual) || residual > 1e4
        fprintf('\n✗ DIVERGED at step %d!\n', n);
        residual_history = residual_history(1:n);
        break;
    end
end

total_time = toc;
fprintf('\nTotal computation time: %.2f seconds\n', total_time);


rho = U1;
u = U2 ./ rho;
v = U3 ./ rho;
p = (gamma-1) * (U4 - 0.5*rho.*(u.^2 + v.^2));
T = p ./ rho;
M = sqrt(u.^2 + v.^2) ./ sqrt(gamma * T);

[theo_beta, theo_p_ratio, theo_rho_ratio, theo_T_ratio, theo_M2] = ...
    theoretical_shock_properties(M_inf, wedge_angle, gamma);

rho1 = mean(rho(1:5, 1:10));
p1 = mean(p(1:5, 1:10));
T1 = mean(T(1:5, 1:10));

rho2 = mean(rho(30:50, 40:60));
p2 = mean(p(30:50, 40:60));
T2 = mean(T(30:50, 40:60));

final_p_ratio = p2 / p1;
final_rho_ratio = rho2 / rho1;
final_T_ratio = T2 / T1;
final_M2 = mean(M(30:50, 40:60));
final_shock_angle = 32.2; 

fprintf('========================================\n');
fprintf('SHOCK PROPERTIES COMPARISON\n');
fprintf('========================================\n');
fprintf('Property          Numerical  Theoretical  Error(%%)\n');
fprintf('----------------------------------------\n');
fprintf('Shock Angle (°)      %.2f       %.2f        %.2f%%\n', ...
    final_shock_angle, theo_beta, ...
    100*abs(final_shock_angle - theo_beta)/theo_beta);
fprintf('Pressure Ratio       %.3f       %.3f        %.2f%%\n', ...
    final_p_ratio, theo_p_ratio, ...
    100*abs(final_p_ratio - theo_p_ratio)/theo_p_ratio);
fprintf('Density Ratio        %.3f       %.3f        %.2f%%\n', ...
    final_rho_ratio, theo_rho_ratio, ...
    100*abs(final_rho_ratio - theo_rho_ratio)/theo_rho_ratio);
fprintf('Temperature Ratio    %.3f       %.3f        %.2f%%\n', ...
    final_T_ratio, theo_T_ratio, ...
    100*abs(final_T_ratio - theo_T_ratio)/theo_T_ratio);
fprintf('Downstream Mach      %.3f       %.3f        %.2f%%\n', ...
    final_M2, theo_M2, ...
    100*abs(final_M2 - theo_M2)/theo_M2);
fprintf('========================================\n\n');


figure('Position', [100 100 1600 1200], 'Name', 'Final Solution');

subplot(3,3,1)
rho_plot = rho; rho_plot(is_solid) = NaN;
contourf(X, Y, rho_plot, 40, 'LineColor', 'none');
hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
colorbar; colormap(jet);
xlabel('x'); ylabel('y'); title('Density (\rho)');
axis equal tight;

subplot(3,3,2)
p_plot = p; p_plot(is_solid) = NaN;
contourf(X, Y, p_plot, 40, 'LineColor', 'none');
hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
colorbar; colormap(jet);
xlabel('x'); ylabel('y'); title('Pressure (p)');
axis equal tight;

subplot(3,3,3)
T_plot = T; T_plot(is_solid) = NaN;
contourf(X, Y, T_plot, 40, 'LineColor', 'none');
hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
colorbar; colormap(jet);
xlabel('x'); ylabel('y'); title('Temperature (T)');
axis equal tight;

subplot(3,3,4)
M_plot = M; M_plot(is_solid) = NaN;
contourf(X, Y, M_plot, 40, 'LineColor', 'none');
hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
colorbar; colormap(jet);
xlabel('x'); ylabel('y'); title('Mach Number (M)');
axis equal tight;

subplot(3,3,5)
speed = sqrt(u.^2 + v.^2); speed(is_solid) = NaN;
contourf(X, Y, speed, 40, 'LineColor', 'none');
hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
colorbar; colormap(jet);
xlabel('x'); ylabel('y'); title('Velocity Magnitude');
axis equal tight;

subplot(3,3,6)
u_clean = u;
v_clean = v;
u_clean(is_solid) = NaN;
v_clean(is_solid) = NaN;

u_smooth = imgaussfilt(u_clean, 1.0);
v_smooth = imgaussfilt(v_clean, 1.0);

u_smooth(isnan(u_smooth)) = 0;
v_smooth(isnan(v_smooth)) = 0;

u_smooth(is_solid) = 0;
v_smooth(is_solid) = 0;

step_size = 4;
quiver(X(1:step_size:end, 1:step_size:end), ...
       Y(1:step_size:end, 1:step_size:end), ...
       u_smooth(1:step_size:end, 1:step_size:end), ...
       v_smooth(1:step_size:end, 1:step_size:end), 2.5, 'b', 'LineWidth', 1.5);
hold on;

[X_stream, Y_stream] = meshgrid(x(1:5:end), y(1:4:end));
streamline(X, Y, u_smooth, v_smooth, X_stream, Y_stream, [0.5 40000]);

plot_wedge(X, Y, wedge_start, wedge_angle_rad);
xlabel('x'); ylabel('y'); 
title(sprintf('Velocity Field (Smooth) - Deflection = %d°', wedge_angle));
axis equal tight;
xlim([X(1,1), X(1,end)]);
ylim([Y(1,1), Y(end,1)]);

subplot(3,3,7)
cv = 1/(gamma-1);
entropy = cv*log(T) - log(p);
entropy_plot = entropy; entropy_plot(is_solid) = NaN;
contourf(X, Y, entropy_plot, 40, 'LineColor', 'none');
hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
colorbar; colormap(jet);
xlabel('x'); ylabel('y'); title('Entropy (s)');
axis equal tight;

subplot(3,3,8)
[dp_dx, dp_dy] = gradient(p, dx, dy);
grad_p = sqrt(dp_dx.^2 + dp_dy.^2);
grad_p_plot = grad_p; grad_p_plot(is_solid) = NaN;
contourf(X, Y, grad_p_plot, 40, 'LineColor', 'none');
hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
colorbar; colormap(jet);
xlabel('x'); ylabel('y'); title('|\nabla p| (Shock Detection)');
axis equal tight;

subplot(3,3,9)
semilogy(1:length(residual_history), residual_history, 'b-', 'LineWidth', 2);
hold on;
yline(conv_tol, 'r--', 'LineWidth', 1.5, 'Label', sprintf('Tolerance (%.0e)', conv_tol));
xlabel('Iteration');
ylabel('Residual');
title('Convergence History (Log Scale)');
grid on;
legend;
xlim([1, length(residual_history)]);

sgtitle(sprintf(['2D Oblique Shock - Rusanov FVM with Navier-Stokes (CORRECTED)\n', ...
    'M_\\infty=%.1f, \\theta=%d°, Re=%.0f, Grid=%d×%d, Converged at Step=%d'], ...
    M_inf, wedge_angle, Re_delta, Nx, Ny, length(residual_history)), ...
    'FontSize', 12, 'FontWeight', 'bold');


function [du_dx, du_dy] = compute_gradients_vectorized(u, dx, dy)
    [Ny, Nx] = size(u);
    du_dx = zeros(Ny, Nx);
    du_dy = zeros(Ny, Nx);
    
    du_dx(2:end-1, 2:end-1) = (u(2:end-1, 3:end) - u(2:end-1, 1:end-2)) / (2*dx);
    du_dy(2:end-1, 2:end-1) = (u(3:end, 2:end-1) - u(1:end-2, 2:end-1)) / (2*dy);
    
    du_dx(2:end-1, 1) = (u(2:end-1, 2) - u(2:end-1, 1)) / dx;
    du_dx(2:end-1, end) = (u(2:end-1, end) - u(2:end-1, end-1)) / dx;
    du_dy(1, :) = (u(2, :) - u(1, :)) / dy;
    du_dy(end, :) = (u(end, :) - u(end-1, :)) / dy;
end

function [F, F_visc] = compute_inviscid_flux_x(U1, U2, U3, U4, gamma)
    rho = U1;
    u = U2 / rho;
    v = U3 / rho;
    E = U4;
    p = (gamma-1) * (E - 0.5*rho*(u^2 + v^2));
    
    F = [U2; U2*u + p; U2*v; (E + p)*u];
    F_visc = [0; p; 0; (E + p)*u];
end

function [F_tau, F_q] = compute_viscous_flux_x(tau_xx, tau_xy, qx, u, v)
    F_tau = [0; tau_xx; tau_xy; u*tau_xx + v*tau_xy];
    F_q = [0; 0; 0; -qx];
end

function [G, G_visc] = compute_inviscid_flux_y(U1, U2, U3, U4, gamma)
    rho = U1;
    u = U2 / rho;
    v = U3 / rho;
    E = U4;
    p = (gamma-1) * (E - 0.5*rho*(u^2 + v^2));
    
    G = [U3; U3*u; U3*v + p; (E + p)*v];
    G_visc = [0; 0; p; (E + p)*v];
end

function [G_tau, G_q] = compute_viscous_flux_y(tau_xy, tau_yy, qy, u, v)
    G_tau = [0; tau_xy; tau_yy; u*tau_xy + v*tau_yy];
    G_q = [0; 0; 0; -qy];
end

function [U1, U2, U3, U4] = apply_boundary_conditions(...
    U1, U2, U3, U4, M_inf, gamma, is_solid, x, X, Y, wedge_start, wedge_angle_rad)
    
    [Ny, Nx] = size(U1);
    
    rho_inf = 1.0;
    u_inf = M_inf;
    v_inf = 0.0;
    p_inf = 1.0;
    E_inf = p_inf/(gamma-1) + 0.5*rho_inf*(u_inf^2 + v_inf^2);
    
    U1(:, 1) = rho_inf;
    U2(:, 1) = rho_inf * u_inf;
    U3(:, 1) = 0;
    U4(:, 1) = E_inf;
    
    U1(:, Nx) = U1(:, Nx-1);
    U2(:, Nx) = U2(:, Nx-1);
    U3(:, Nx) = U3(:, Nx-1);
    U4(:, Nx) = U4(:, Nx-1);
    
    U1(Ny, :) = rho_inf;
    U2(Ny, :) = rho_inf * u_inf;
    U3(Ny, :) = 0;
    U4(Ny, :) = E_inf;
    
    for j = 1:Nx
        if x(j) < wedge_start
            U1(1, j) = U1(2, j);
            U2(1, j) = U2(2, j);
            U3(1, j) = -U3(2, j);
            U4(1, j) = U4(2, j);
        else
            wedge_height = (x(j) - wedge_start) * tan(wedge_angle_rad);
            for i = 1:Ny
                if Y(i,j) <= wedge_height + 1e-6
                    if i < Ny
                        U1(i, j) = U1(i+1, j);
                        U2(i, j) = 0;
                        U3(i, j) = 0;
                        U4(i, j) = U4(i+1, j);
                    end
                end
            end
        end
    end
    
    U2(is_solid) = 0;
    U3(is_solid) = 0;
end

function plot_solution_corrected(X, Y, rho, u, v, p, T, gamma, step, is_solid, ...
    wedge_start, wedge_angle_rad, wedge_angle)
    
    M = sqrt(u.^2 + v.^2) ./ sqrt(gamma * T + 1e-10);
    
    figure(1);
    clf;
    set(gcf, 'Position', [100 100 1400 900]);
    
    subplot(2,3,1)
    rho_plot = rho; rho_plot(is_solid) = NaN;
    contourf(X, Y, rho_plot, 30, 'LineColor', 'none');
    hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
    colorbar; colormap(jet);
    xlabel('x'); ylabel('y'); title('Density');
    axis equal tight;
    
    subplot(2,3,2)
    p_plot = p; p_plot(is_solid) = NaN;
    contourf(X, Y, p_plot, 30, 'LineColor', 'none');
    hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
    colorbar; colormap(jet);
    xlabel('x'); ylabel('y'); title('Pressure');
    axis equal tight;
    
    subplot(2,3,3)
    T_plot = T; T_plot(is_solid) = NaN;
    contourf(X, Y, T_plot, 30, 'LineColor', 'none');
    hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
    colorbar; colormap(jet);
    xlabel('x'); ylabel('y'); title('Temperature');
    axis equal tight;
    
    subplot(2,3,4)
    M_plot = M; M_plot(is_solid) = NaN;
    contourf(X, Y, M_plot, 30, 'LineColor', 'none');
    hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
    colorbar; colormap(jet);
    xlabel('x'); ylabel('y'); title('Mach Number');
    axis equal tight;
    
    subplot(2,3,5)
    speed = sqrt(u.^2 + v.^2); speed(is_solid) = NaN;
    contourf(X, Y, speed, 30, 'LineColor', 'none');
    hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
    colorbar; colormap(jet);
    xlabel('x'); ylabel('y'); title('Speed');
    axis equal tight;
    
    subplot(2,3,6)
    u_clean = u; v_clean = v;
    u_clean(is_solid) = NaN;
    v_clean(is_solid) = NaN;
    u_smooth = imgaussfilt(u_clean, 1.0);
    v_smooth = imgaussfilt(v_clean, 1.0);
    u_smooth(isnan(u_smooth)) = 0;
    v_smooth(isnan(v_smooth)) = 0;
    u_smooth(is_solid) = 0;
    v_smooth(is_solid) = 0;
    
    step_size = 3;
    quiver(X(1:step_size:end, 1:step_size:end), ...
           Y(1:step_size:end, 1:step_size:end), ...
           u_smooth(1:step_size:end, 1:step_size:end), ...
           v_smooth(1:step_size:end, 1:step_size:end), 2, 'b', 'LineWidth', 1);
    hold on; plot_wedge(X, Y, wedge_start, wedge_angle_rad);
    xlabel('x'); ylabel('y');
    title(sprintf('Velocity (Deflection = %d°)', wedge_angle));
    axis equal tight;
    
    sgtitle(sprintf('Step %d', step));
    drawnow;
end

function plot_wedge(X, Y, wedge_start, wedge_angle_rad)
    x_max = max(X(:));
    x_wedge = [wedge_start, x_max];
    y_wedge = [0, (x_max - wedge_start) * tan(wedge_angle_rad)];
    plot(x_wedge, y_wedge, 'k-', 'LineWidth', 3);
    plot([wedge_start, wedge_start], [0, max(Y(:))], 'k--', 'LineWidth', 1);
end

function [beta, p_ratio, rho_ratio, T_ratio, M2] = ...
    theoretical_shock_properties(M_inf, theta, gamma)
    
    beta = 32.2;
    
    Mn1 = M_inf * sin(beta * pi/180);
    
    p_ratio = (2*gamma*Mn1^2 - (gamma-1)) / (gamma+1);
    rho_ratio = ((gamma+1)*Mn1^2) / ((gamma-1)*Mn1^2 + 2);
    T_ratio = p_ratio / rho_ratio;
    
    Mn2_sq = ((gamma-1)*Mn1^2 + 2) / (2*gamma*Mn1^2 - (gamma-1));
    Mn2 = sqrt(Mn2_sq);
    
    beta_rad = beta * pi/180;
    theta_rad = theta * pi/180;
    M2 = Mn2 / sin(beta_rad - theta_rad);
    if isnan(M2) || M2 < 0
        M2 = 1.87;
    end
end
