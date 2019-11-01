addpath ServiceRoutines
% Driver script for solving the 1D(spherically symmetric reduction) Generalized harmonic Einstein equations
Globals1D;

% Order of polymomials used for approximation 
N = 2;

% Generate simple mesh
[Nv, VX, K, EToV] = MeshGen1D(1.5, 10, 50);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
%init_min;
init_schwarzschild_kerr_schild;
%init_PG;
%init_test;
time = 0;                                                                                                
FinalTime = 10000.0; 
time_seq = [];
rhs_psi_seq = [];
rhs_Pi_psi_seq = [];
rhs_Phi_psi_seq = [];
Cr_psi_seq = [];
                                                                                                         
% Runge-Kutta residual storage                                                                           
res_psi = zeros(Np,K);                                                                                      
res_Pi_psi = zeros(Np,K);                                                                                      
res_Phi_psi = zeros(Np,K);                                                                                      
                                                                                                         
% compute time step size                                                                                 
xmin = min(abs(x(1,:)-x(2,:)));                                                                          
dt = 0.9/(2*N+1)*xmin;                                                            
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;                                                      

b = sqrt(gamma11(vmapO))
for tstep=1:Nsteps                                                                                   
    for INTRK = 1:5                                                                            
        timelocal = time + rk4c(INTRK)*dt;                                                           

        [rhs_psi, rhs_Pi_psi, rhs_Phi_psi] = compute_RHS;
        rhs_Pi_psi(vmapO) = paragamma2/2*rhs_psi(vmapO) + 1/2*rhs_Pi_psi(vmapO) ...
                + b/2*rhs_Phi_psi(vmapO);
        rhs_Phi_psi(vmapO) = -paragamma2/2/b*rhs_psi(vmapO) ...
                + 1/2/b*rhs_Pi_psi(vmapO) + 1/2*rhs_Phi_psi(vmapO);

        res_psi = rk4a(INTRK)*res_psi + dt*rhs_psi;                                                           
        res_Pi_psi = rk4a(INTRK)*res_Pi_psi + dt*rhs_Pi_psi;                                                           
        res_Phi_psi = rk4a(INTRK)*res_Phi_psi + dt*rhs_Phi_psi;                                                           

        psi = psi+rk4b(INTRK)*res_psi;                                                                      
        Pi_psi = Pi_psi+rk4b(INTRK)*res_Pi_psi;                                                                      
        Phi_psi = Phi_psi+rk4b(INTRK)*res_Phi_psi;                                                                      
    end;                                                                                             

    % limiter
    %g00 = SlopeLimitN(g00);
    %g01 = SlopeLimitN(g01);
    %g11 = SlopeLimitN(g11);

    %filter
    %F = Filter1D(N, 0, 50);
    %g00 = F*g00;
    %g01 = F*g01;
    %Phi11 = F*Phi11;

    % Increment time                                                                                 
    time = time+dt;                                                                                  
    time_seq = [time_seq, time];
    rhs_psi_seq = [rhs_psi_seq, max(max(abs(rhs_psi)))];
    rhs_Pi_psi_seq = [rhs_Pi_psi_seq, max(max(abs(rhs_Pi_psi)))];
    rhs_Phi_psi_seq = [rhs_Phi_psi_seq, max(max(abs(rhs_Phi_psi)))];

    Cr_psi_seq = [Cr_psi_seq, max(max(abs(Cr_psi)))];
    if (mod(tstep, 20) == 0)
        figure(1); plot(x, psi); title(['psi, t = ', num2str(time)]); drawnow; pause(.1);
        %figure(2); plot(x, Pi_psi); title(['Pi\_psi, t = ', num2str(time)]); drawnow; pause(.1);
        %figure(3); plot(x, Phi_psi); title(['Phi\_psi, t = ', num2str(time)]); drawnow; pause(.1);

        figure(4); plot(x, Cr_psi); title(['Cr\_psi, t = ', num2str(time)]); drawnow; pause(.1);

        %figure(6); semilogy(time_seq, Cr_psi_seq); title(['max of Cr\_psi with time']); drawnow; pause(.1);
        %figure(7); semilogy(time_seq, rhs_psi_seq); title(['max of rhs\_psi with time']); drawnow; pause(.1);
        %figure(8); semilogy(time_seq, rhs_Pi_psi_seq); title(['max of rhs\_Pi\_psi with time']); drawnow; pause(.1);
        %figure(9); semilogy(time_seq, rhs_Phi_psi_seq); title(['max of rhs\_Phi\_psi with time']); drawnow; pause(.1);

    end;
end;  

