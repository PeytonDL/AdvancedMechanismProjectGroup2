%%  Genetic Algorithm
clear; close; 
dtor = (pi / 180); %Degree-to-radian conversion factor
rtod = (180 / pi); %Radian-to-degree conversion factor

%% Problem Definition

% % X-Y coordinates of Prescribed Positions
pps(1) = 0 + 0.0*1i;   % PP 1 coords
pps(2) =  1 + 1i;  % PP 2 coords 
pps(3) =  2.5 + 1i;   % PP 3 coords 
pps(4) = 4 + 0*1i;   % PP 4 coords 
betaD1(1) = 0;
betaD2(1) = 0;

betaD1(2) = -25*dtor; %Initial choice for Beta 2, converted to radians
betaD1(3) = -55*dtor;

betaD2(2) = 70*dtor; %Second dyad choice of Beta 2, converted to radians
betaD2(3) = 35*dtor;

alpha(1) = 0;
alpha(2) = 5*dtor;
alpha(3) = 50*dtor;% Alpha values are shared by both dyads

% Displacement vector calculation if using PPs
delta(2) = pps(2) - pps(1); % PP displacement from PP1 to PP2
delta(3) = pps(3) - pps(1); % PP displacement from PP1 to PP3
delta(4) = pps(4) - pps(1); % PP displacement from PP1 to PP2

%% Genetic/Evolutionary Algorithm Parameters
populationSize = 50;
numGenerations = 50;
mutationRate = 0.05;
fitscore = 0;
generations = zeros(numGenerations,12);

%% Synthesize Initial Solution
A=[exp(1i*betaD1(2))-1 exp(1i*alpha(2))-1;
   exp(1i*betaD1(3))-1 exp(1i*alpha(3))-1];

deltavect=[delta(2);
           delta(3)];

Sol_Dyad1=A\deltavect;
W1D1=Sol_Dyad1(1);
Z1D1=Sol_Dyad1(2);
OAD1=pps(1)-Z1D1-W1D1;
%% Second Dyad
A=[exp(1i*betaD2(2))-1 exp(1i*alpha(2))-1;
   exp(1i*betaD2(3))-1 exp(1i*alpha(3))-1];

deltavect=[delta(2);
           delta(3)];

Sol_Dyad1=A\deltavect;
W1D2=Sol_Dyad1(1);
Z1D2=Sol_Dyad1(2);
OAD2=pps(1)-Z1D2-W1D2;

generations(1,:) = [W1D1, Z1D1, OAD1, betaD1(2), betaD1(3), W1D2, Z1D2, OAD2, betaD2(2), betaD2(3), fitscore, 1]; %Dyad 1 data, Dyad 2 data, fitness score, generation number


text = sprintf('Initial Dyads');
figure('Position', [10 50 800 600], 'Name', text,'NumberTitle', 'off'); 
hold on;
axis equal;
grid on;
title('Initial Solution');
h=axis;
scale=(h(2)-h(1));

for j = 1: 3
  
     %Plot WA
     quiver(real(OAD1),imag(OAD1),real(W1D1*(exp(betaD1(j)*1i))),   imag(W1D1 *(exp(betaD1(j)*1i))),1, 'Color', 'c','LineWidth',2);
     %Plot ZA
     quiver(real(OAD1+W1D1*(exp(betaD1(j)*1i))), imag(OAD1+W1D1*(exp(betaD1(j)*1i))),real(Z1D1*(exp(alpha(j)*1i))),  imag(Z1D1*(exp(alpha(j)*1i))),1, 'Color', 'b','LineWidth',2);    

    
     %Plot WA
     quiver(real(OAD2),imag(OAD2),real(W1D2*(exp(betaD2(j)*1i))),   imag(W1D2 *(exp(betaD2(j)*1i))),1, 'Color', 'r','LineWidth',2);
     %Plot ZA
     quiver(real(OAD2+W1D2*(exp(betaD2(j)*1i))), imag(OAD2+W1D2*(exp(betaD2(j)*1i))),real(Z1D2*(exp(alpha(j)*1i))),  imag(Z1D2*(exp(alpha(j)*1i))),1, 'Color', 'g','LineWidth',2);   

     plot([real(OAD2+W1D2*(exp(betaD2(j)*1i))), real(OAD1+W1D1*(exp(betaD1(j)*1i)))],[imag(OAD2+W1D2*(exp(betaD2(j)*1i))),imag(OAD1+W1D1*(exp(betaD1(j)*1i)))],'k','LineWidth',2);
% 
end
quiver(real(pps(1)),imag(pps(1)),1,0,scale,'ok');
quiver(real(pps(1))+real(delta(2)),imag(pps(1))+imag(delta(2)),real(exp(gamma(2)*1i)),imag(exp(gamma(2)*1i)),scale,'ok'); 
quiver(real(pps(1))+real(delta(3)),imag(pps(1))+imag(delta(3)),real(exp(gamma(3)*1i)),imag(exp(gamma(3)*1i)),scale,'ok');
quiver(real(pps(1))+real(delta(4)),imag(pps(1))+imag(delta(4)),real(exp(gamma(4)*1i)),imag(exp(gamma(4)*1i)),scale,'ok');

axis equal;
hold off;

%% Run Calculations
% Initialization
population = initializePopulation(generations(1,:),populationSize, mutationRate, alpha, delta, pps);

for generation = 1:numGenerations
    % Evaluation
    fitness = evaluatePopulation(population);
    
    % Selection
    selectedParent = selection(fitness);

    % Replace population with new generation
    population = initializePopulation(selectedParent,populationSize, mutationRate, alpha, delta, pps);
    % Display best individual in current generation
%     [~, bestIndex] = max(fitness);
%     bestIndividual = population(bestIndex,:);
%     disp(['Generation ', num2str(generation), ', Best Fitness: ', num2str(fitness(bestIndex))]);
end
%% Plot Final Solution
%[W1D1, Z1D1, OAD1, betaD1(2), betaD1(3), W1D2, Z1D2, OAD2, betaD2(2), betaD2(3), fitscore, 1]; %Dyad 1 data, Dyad 2 data, fitness score, generation number

Fitscore = selectedParent(11)*rtod
W1D1 = selectedParent(1);
Z1D1 = selectedParent(2);
OAD1 = selectedParent(3);
betaD1(2) = selectedParent(4);
betaD1(3) = selectedParent(5);
W1D2 = selectedParent(6);
Z1D2 = selectedParent(7);
OAD2 = selectedParent(8);
betaD2(2) = selectedParent(9);
betaD2(3) = selectedParent(10);

text = sprintf('Solution Dyads');
figure('Position', [10 50 800 600], 'Name', text,'NumberTitle', 'off'); 
hold on;
axis equal;
grid on;
title('Solution Optimized for Transmission Angle');
h=axis;
scale=(h(2)-h(1))/10;

for j = 1: 3
  
     %Plot WA
     quiver(real(OAD1),imag(OAD1),real(W1D1*(exp(betaD1(j)*1i))),   imag(W1D1 *(exp(betaD1(j)*1i))),1, 'Color', 'c','LineWidth',2);
     %Plot ZA
     quiver(real(OAD1+W1D1*(exp(betaD1(j)*1i))), imag(OAD1+W1D1*(exp(betaD1(j)*1i))),real(Z1D1*(exp(alpha(j)*1i))),  imag(Z1D1*(exp(alpha(j)*1i))),1, 'Color', 'b','LineWidth',2);    

    
     %Plot WA
     quiver(real(OAD2),imag(OAD2),real(W1D2*(exp(betaD2(j)*1i))),   imag(W1D2 *(exp(betaD2(j)*1i))),1, 'Color', 'r','LineWidth',2);
     %Plot ZA
     quiver(real(OAD2+W1D2*(exp(betaD2(j)*1i))), imag(OAD2+W1D2*(exp(betaD2(j)*1i))),real(Z1D2*(exp(alpha(j)*1i))),  imag(Z1D2*(exp(alpha(j)*1i))),1, 'Color', 'g','LineWidth',2);   

     plot([real(OAD2+W1D2*(exp(betaD2(j)*1i))), real(OAD1+W1D1*(exp(betaD1(j)*1i)))],[imag(OAD2+W1D2*(exp(betaD2(j)*1i))),imag(OAD1+W1D1*(exp(betaD1(j)*1i)))],'k','LineWidth',2);
end
quiver(real(pps(1)),imag(pps(1)),1,0,scale,'ok');
quiver(real(pps(1))+real(delta(2)),imag(pps(1))+imag(delta(2)),real(exp(gamma(2)*1i)),imag(exp(gamma(2)*1i)),scale,'ok'); 
quiver(real(pps(1))+real(delta(3)),imag(pps(1))+imag(delta(3)),real(exp(gamma(3)*1i)),imag(exp(gamma(3)*1i)),scale,'ok');

axis equal;
hold off;



%% Functions Definitions




function population = initializePopulation(currentgen, populationSize, mutationRate, alpha, delta, pps)
    % Initialize population with random solutions
    %Define limits to mutation
    % Define mutations: 
    fitscore = 0;
    
    population = zeros(populationSize,12);
    for i = 1:1:populationSize
        
        if(rand(1,1)<mutationRate) %Check to mutate beta 2 of dyad 1
             betaD1(2) = currentgen(4)+(1+1)*rand(1,1) - 1;
        else
            betaD1(2) = currentgen(4);
        end
        if(rand(1,1)<mutationRate) %Check to mutate beta 3 of dyad 1
             betaD1(3) = currentgen(5)+(1+1)*rand(1,1) - 1;
        else
            betaD1(3) = currentgen(5);
        end
        if(rand(1,1)<mutationRate) %Check to mutate beta 2 of dyad 2
             betaD2(2) = currentgen(9)+(1+1)*rand(1,1) - 1;
        else
            betaD2(2) = currentgen(9);
        end
        if(rand(1,1)<mutationRate) %Check to mutate beta 3 of dyad 2
             betaD2(3) = currentgen(10)+(1+1)*rand(1,1) - 1;
        else
            betaD2(3) = currentgen(10);
        end
        
        A=[exp(1i*betaD1(2))-1 exp(1i*alpha(2))-1;
           exp(1i*betaD1(3))-1 exp(1i*alpha(3))-1];

        deltavect=[delta(2);
                   delta(3)];

        Sol_Dyad1=A\deltavect; %Calculate W and Z
        W1D1=Sol_Dyad1(1);
        Z1D1=Sol_Dyad1(2);
        OAD1=pps(1)-Z1D1-W1D1;

        A=[exp(1i*betaD2(2))-1 exp(1i*alpha(2))-1;
           exp(1i*betaD2(3))-1 exp(1i*alpha(3))-1];

        deltavect=[delta(2);
                   delta(3)];

        Sol_Dyad1=A\deltavect; %Calculate W and Z
        W1D2=Sol_Dyad1(1);
        Z1D2=Sol_Dyad1(2);
        OAD2=pps(1)-Z1D2-W1D2;
        population(i,:) = [W1D1, Z1D1, OAD1, betaD1(2), betaD1(3), W1D2, Z1D2, OAD2, betaD2(2), betaD2(3), fitscore, currentgen(12)];
    end 

end



function fitness = evaluatePopulation(population)
    % Evaluate fitness of each individual in the population
    
    for i = 1:1:length(population)
        L1 = abs(population(i,3) + population(i,1) + population(i,2) - population(i,7) - population(i,6));
        L2 = abs(population(i,1));
        L3 = abs(population(i,2) - population(i,7));
        L4 = abs(population(i,6));
        theta2 = angle(population(i,1));
        theta3 = angle(population(i,2)-population(i,7));
        transmission_ang = atan2(sqrt(1-(((L2^2 + L3^2 - L1^2 - L4^2) - 2*L2*L3*cos(theta2 - theta3))/(2*L2*L3)))^2, ((L2^2 + L3^2 - L1^2 - L4^2) - 2*L2*L3*cos(theta2 - theta3))/(2*L2*L3)); 
        population(i,11) = transmission_ang;
    end
    fitness = population;
end



function selectParent = selection(population)
    % Select parent based on fitness
    fitmax = [1,0];
    for i=1:1:length(population)    
        if(population(i,11)>fitmax(2))
            fitmax(1)=i;
            fitmax(2)=population(i,11);
        end
    end
    selectParent = population(fitmax(1),:);
    %score = population(fitmax(1),11)
end

% function offspring = crossover(selectedParents)
%     % Single-point crossover
%     crossoverPoint = randi([1, size(selectedParents, 2) - 1]);
%     offspring = [selectedParents(1,1:crossoverPoint), selectedParents(2,crossoverPoint+1:end);
%                  selectedParents(2,1:crossoverPoint), selectedParents(1,crossoverPoint+1:end)];
% end
% 
% function mutatedOffspring = mutation(offspring, mutationRate)
%     % Bit-wise mutation
%     mutationMask = rand(size(offspring)) < mutationRate;
%     mutatedOffspring = xor(offspring, mutationMask);
% end

