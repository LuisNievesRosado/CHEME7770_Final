function out =  MC(Ksub,eta)

%=========================================================================%
% Kinetic Monte Carlo Code for the molecular clutch mechanism with a      
% viscoelastic substrate. Based on method by Chan and Odde                
% DOI:10.1126/science.1163595                                             
%=========================================================================%

%=========================================================================%
% Parameter Initialization
%=========================================================================%
Fm = -150;      % Total myosin pulling force [pN]
nuu = -120;     % Unloaded motor velocity [nm/s]
nc = 75;        % Total Number of clutches []
kon = 1;        % Clutch attachment rate constant [s^-1]
koff0 = 0.1;    % Clutch detachment rate constant [s^-1]
Fb = -2;        % Charactheric bond rupture force [pN]
Kclutch = 5;    % Clutch spring constant [pN/nm]
%Ksub = 5;      % Substrate Spring Constant [pN/nm]
%eta = 0.05;        % Substrate viscosity [pN s/nm]
%=========================================================================%


%=========================================================================%
% System Variable Initialization
%=========================================================================%
dt = 0.005;                % Timestep. 5ms
N = 100000;                 % Number of timeteps

xclutch = zeros(1,nc);      % Location of clutches
Fclutch = zeros(1,nc);      % Forces of clutches
clutchstatus = zeros(1,nc); % Status of clutch. 0 =  detached 1 = attached
koff = zeros(1,nc);         % Strain dependent detachment rate

xsub = zeros(1,N);          % Location of substrate
xdot = 0;                   % Strain rate of substrate

nu = zeros(1,N);            % Actin speed

Na = sum(clutchstatus);     % Number of attached clutches
Nd = nc - Na;               % Number of detached clutches

%=========================================================================%
% Simulation
%=========================================================================%

for i = 1:N
    % Compute koff for all clutches
    for j = 1:nc
        Fclutch(j) = Kclutch*(xclutch(j)-xsub(i));
        koff(j) = koff0*exp(Fclutch(j)/Fb);
    end
    
    % Compute Rtot
    Rtot = Nd*kon + sum(koff.*clutchstatus);
    
    % Determine attachments and deattachments
    for j = 1:nc
       if clutchstatus(j) == 1           % Clutch is attached, can deattach
           lim = koff(j)/Rtot;
           check = rand(1);
           if check < lim                % Clutch deataches
               clutchstatus(j) = 0;
               %xclutch(j) = 0;
           end
       elseif clutchstatus(j) == 0       % Clutch is detached, can attach
           lim = kon/Rtot;
           check = rand(1);
           if check < lim                % Clutch attaches
               clutchstatus(j) = 1;
           end
       end
    end
    
    % Determine clutch counts
    Na = sum(clutchstatus);             % Number of attached clutches
    Nd = nc - Na;                       % Number of detached clutches
    
    % Change substrate position
    xold = xsub(i);
    xsub(i) = (Kclutch*sum(xclutch.*clutchstatus))/(Ksub + eta*xdot + Kclutch*Na);
    
    % Update xdot
    xdot = (xsub(i)-xold)/dt;
    
    % Compute actin velocity
    nu(i) = nuu*(1-(Ksub*xsub(i)+eta*xdot)/Fm);
    
    % Update attached clutch location
    for j = 1:nc
        if clutchstatus(j) == 1         % Clutch is attached
            xclutch(j) = xclutch(j) + nu(i)*dt;
        end
    end
end

%=========================================================================%
% Data Extraction
%=========================================================================%
% Plotting
% time = 0:dt:dt*(N-1);
% figure;
% plot(time,nu)
% title('Speed')
% figure;
% plot(time,xsub)
% title('Location')

out = mean(nu(N-20000:N));


