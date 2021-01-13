%% PSO With Variable Accelerator Coefficients

% Author: Dr. Davide Conte
% Contact Info: davide.conte90@gmail.com


%% Clearing memory, screen, etc.
close all; clear all; clc;

global P J JBest PBest GG N_particles N_elements V BLv BUv BLp BUp
global N_iterations GBest
% useful if making any of the following parts into functions

%-------------------------------------------------------------------------%
% choose one of the examples (1 or 2)
Example = 1;
% 1 --> f(x,y,z) = (x-1)^2 + (y-2)^2 + (z-3)^2
% 2 --> f(x,y) = 3x^2 + 2xy + y^2 - 4x + 5y
%-------------------------------------------------------------------------%

%% Basic PSO initialization

N_particles = 30;
% [ ] number of particles to initiate

if Example == 1
    N_elements = 3;
    % [ ] number of elements each particle has (i.e. parameters)
elseif Example == 2
    N_elements = 2;
end

N_iterations = 500;
% [ ] number of maximum iterations before PSO is stopped


%% Particle initialization, particle "position" and "velocity"

BLp = -20.*ones(N_elements); BUp = 20.*ones(N_elements);
% [ ] set lower and upper bounds on unknowns (particle elements)

% create particles using uniform distribution according to BLp and BUp
for k = 1:1:N_elements
    
    P = BLp(k).*ones(N_particles,k)+ ...
        (BUp(k)-BLp(k)).*rand(N_particles,k); 
    
    % could also use the built-in function "unifrnd"

end

PBest = zeros(N_particles,N_elements);
% [ ] initialize the "best" for each particle to all zeros (arbitrary)

GGstar = zeros(1,N_iterations);
% [ ] initialize the vector of best J's for all iterations

V = zeros(N_particles, N_elements);
% [ ] initialize particle velocity for each particle

BUv = BUp - BLp; BLv = -BUv;
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

for j = 1:N_iterations
    
    %---------------------------------------------------------------------%    
    % Evaluate cost function J for each particle in current iteration
    for i = 1:N_particles
        if Example == 1
            J(i) = (P(i,1)-1)^2 + (P(i,2)-2)^2 + (P(i,3)-3)^2;
            % f(x,y,z) = (x-1)^2 + (y-2)^2 + (z-3)^2 (example 1)
        elseif Example == 2
            J(i) = 3*P(i,1)^2 + 2*P(i,1)*P(i,2) + P(i,2)^2 - ...
                4*P(i,1) + 5*P(i,2);
            % f(x,y) = 3x^2 + 2xy + y^2 - 4x + 5y (example 2)
        end
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
        clc; fprintf('Progress: %2.0f %%\n',j/N_iterations*100)
    end
    
end

PSOtime = toc;
% [s] time it takes to run PSO 

time = datestr(now, 'yyyy_mm_dd');
filename = sprintf('PSO_%s.mat',time);
save( fullfile('', filename) )
% save workspace with today's date


%% Final Results
fprintf('\nThe particle with best J is: \n'); disp(GBest)
fprintf('\nAnd its corresponding best J is: \n'); disp(GG)

% output run time to the screen 
if PSOtime > 3600
    fprintf('\nThe PSO algorithm took %2.2f hours to complete\n',PSOtime)
elseif PSOtime > 60
    fprintf('\nThe PSO algorithm took %2.2f minutes to complete\n',PSOtime)
else
    fprintf('\nThe PSO algorithm took %2.4f seconds to complete\n',PSOtime)
end

%figure of J vs. iteration number (semilog plot if J > 0 for all J's)
if GGstar(end) > 0
    figure; semilogy(1:N_iterations, GGstar)
    title('J vs. Iteration Number')
    xlabel('Iteration Number'); ylabel('J (Cost Function)');
    grid on; set(gca,'Fontsize',20);
else
    figure; plot(1:N_iterations, GGstar)
    title('J vs. Iteration Number')
    xlabel('Iteration Number'); ylabel('J (Cost Function)');
    grid on; set(gca,'Fontsize',20);
end


% script ends here