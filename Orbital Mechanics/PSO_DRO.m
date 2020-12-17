%% PSO With Variable Accelerator Coefficients

% Author: Dr. Davide Conte
% Contact Info: davide.conte90@gmail.com
% Co-Author: Vittorio Baraldi
% Contact: baraldiv@my.erau.edu

%% Clearing memory, screen, etc.
close all; clear all; clc;

global P J JBest PBest GG N_particles N_elements V BLv BUv BLp BUp
global N_iterations GBest lambda a b
% useful if making any of the following parts into functions

%-------------------------------------------------------------------------%
%% Initial data
mu_m=42828.375214;              %mars' gravitational parameter [km^3/s^2]
mu_p=7.113588120963051e-4;      %phobos' gravitational parameter [km^3/s^2]
R=9376;                         %average distance from m1 and m2 [km]
lambda=mu_p/(mu_m+mu_p);        %mass parameter
LU=R;
TU=sqrt(R^3/(mu_m+mu_p));
a=13/LU;                        %semi-major axis of phobos approximating ellipsis
b=11.4/LU;                      %semi-minor axis of phobos approximating ellipsis

%Given DRO for AX=100km
Vy=-0.045620256764708;
T=2.731044880670166e4;

%find DRO for following Ax
Ax=125;
%nondimensional parameters
AX=Ax/LU;
T1=T/TU;
VY=Vy*TU/LU;
%-------------------------------------------------------------------------%
%% Calculating Phobos Approximation Ellipse
xph1=-a/2+(1-lambda);
xph2=a/2+(1-lambda);
t=linspace(0,2*pi,100);
X=a*cos(t);
Y=b*sin(t);
w = atan2(0,xph2-xph1);
xph = (xph1+xph2)/2 + X*cos(w) - Y*sin(w);
yph = X*sin(w) + Y*cos(w);

%% Basic PSO initialization

N_particles = 40;
% [ ] number of particles to initiate

N_elements = 2;

N_iterations = 500;
% [ ] number of maximum iterations before PSO is stopped


%% Particle initialization, particle "position" and "velocity"

BLp(1,1) = 1.5*VY; BUp(1,1) = VY;
BLp(2,1) = T1; BUp(2,1) = 2*pi;
% [ ] set lower and upper bounds on unknowns (particle elements)

% create particles using uniform distribution according to BLp and BUp
for k = 1:N_particles
    for j=1:N_elements
        P(k,j)=unifrnd(BLp(j,1),BUp(j,1));
        % could also use the built-in function "unifrnd"
    end
end

PBest = zeros(N_particles,N_elements);
% [ ] initialize the "best" for each particle to all zeros (arbitrary)

GGstar = zeros(1,N_iterations);
% [ ] initialize the vector of best J's for all iterations

V = zeros(N_particles, N_elements);
% [ ] initialize particle velocity for each particle

BUv = BUp-BLp; BLv = -BUv;
% [ ] determine velocity bounds from position

c_I = (1 + rand)/2; 
c_C = 1.49445*rand;
c_S = 1.49445*rand;
% [ ] acceleration coefficients (vary with application)


%% Cost function

J = zeros(N_particles);
for kk = 1:N_particles
    JBest(kk) = inf;
end
GG = inf;
% [ ] initialization of cost function (J), best overall particle (JBest),
% and overall best cost function obtained (GG)
% NB: change "inf" to "-inf" if maximizing J instead of minimizing


%% Main PSO computations

tic; 
% keep track of time it takes to run PSO computations (use profiler to see
% where improvements can be made

%initializing r0 vector
r0=[AX+(1-lambda);0];
for j = 1:N_iterations
 
    %---------------------------------------------------------------------%    
    % Evaluate cost function J for each particle in current iteration
    for i = 1:N_particles
        
    %integrating the CR3BP EOMs
    v0=[0;P(i,1)];
    y0=[r0;v0];
    timespan=[0 P(i,2)];
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,yf]=ode113(@orbitfunCR3BP,timespan,y0,options);
    pos=[yf(end,1);yf(end,2)];
    vel=[yf(end,3);yf(end,4)];
    %evaluating cost function
    J(i)= norm(pos-r0)+norm(vel-v0);
    end
    %---------------------------------------------------------------------%
    
    % Evaluate PBest and GBest (best particle and best cost)
    for i = 1:N_particles
        if J(i) < JBest(i)
            PBest(i,:) = P(i,:);
            JBest(i) = J(i);
        end
    end
    
    for i = 1:N_particles
        if J(i) < GG
           GG = J(i);
           GBest = P(i,:);
        end
    end
    
    
    % Update particle velocity 
    for i = 1:N_particles
        
        V(i,:) = c_I*V(i,:) + c_C*(PBest(i,:) - P(i,:)) + ...
            c_S*(GBest - P(i,:));
        % new particle velocity 
        
        % check if boundaries are being respected
        for k = 1:N_elements
            if V(i,k) < BLv(k)
                V(i,k) = BLv(k);
            end
            if V(i,k) > BUv(k)
                V(i,k) = BUv(k);
            end
        end
    end
    
    
	% update particle position for each particle
    for i = 1:N_particles
        
        P(i,:) = P(i,:) + V(i,:);
        % new particle position
        
        % check position boundaries
        for k = 1:N_elements
            if P(i,k) < BLp(k)
                P(i,k) = BLp(k);
                V(i,k) = 0;
            end
            if P(i,k) > BUp(k)
                P(i,k) = BUp(k);
                V(i,k) = 0;
            end
        end
    end
    
    GGstar(j) = GG;
    
    % give progress updates to the screen every 10%
    if mod(j/N_iterations*100,10) == 0
        clc; fprintf('Progress of PSO for a DRO (Ax = %i km): %2.0f %%\n', Ax,j/N_iterations*100)
    end
    
end
PSOtime = toc;
% [s] time it takes to run PSO 

time = datestr(now, 'yyyy_mm_dd');
filename = sprintf('PSO_DRO_%s.mat',time);
save( fullfile('', filename) )
% save workspace with today's date

%% Plot landing trajectory corresponding to best result
figure(1)
hh(1)=plot(yf(:,1)*LU,yf(:,2)*LU,'LineWidth',1.5); hold on;
hh(2)=plot(xph*LU,yph*LU,'Color',[0.5,0.5,0.5],'LineWidth',1.5); hold on;
hh(3)=patch(xph*LU,yph*LU,[0.5,0.5,0.5]); hold on
axis equal

title('Mars-Phobos DRO'); xlabel('x-rotating [km]'); 
ylabel('y-rotating [km]'); hold on; grid on; axis equal; 
set(gca,'FontSize',20);

legend(hh([1 3]),'DRO','Phobos');

%% Final Results
fprintf('\nThe particle with best J is: \n'); disp(GBest)
fprintf('\nAnd its corresponding best J is: \n'); disp(GG)
fprintf('The y component of the velocity for Ax = %i is:\n%1.8f km/s\n',Ax,GBest(1,1)*LU/TU)
fprintf('The orbital period for Ax = %i is:\n%1.8f s\n',Ax,GBest(1,2)*TU)
finalVY=GBest(1,1)*LU/TU;
finalT=GBest(1,2)*TU;
% output run time to the screen 
if PSOtime > 3600
    fprintf('\nThe PSO algorithm took %2.2f hours to complete\n',PSOtime/3600)
elseif PSOtime > 60
    fprintf('\nThe PSO algorithm took %2.2f minutes to complete\n',PSOtime/60)
else
    fprintf('\nThe PSO algorithm took %2.4f seconds to complete\n',PSOtime)
end

%figure of J vs. iteration number (semilog plot if J > 0 for all J's)
if GGstar(end) > 0
    figure; semilogy(1:N_iterations, GGstar)
    xlabel('Iteration Number'); ylabel('J (Cost Function)');
    grid on; set(gca,'Fontsize',30);
else
    figure; plot(1:N_iterations, GGstar)
    xlabel('Iteration Number'); ylabel('J (Cost Function)');
    grid on; set(gca,'Fontsize',30);
end


% script ends here