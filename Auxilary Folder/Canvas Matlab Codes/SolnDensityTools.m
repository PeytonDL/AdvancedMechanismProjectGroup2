%M+K Circle Plotting, Dyad Pair Solution
%
% Authors: Sean Mather
%
% Date:    July 28, 2025


% Create a density plot %

clear; clc; close all;

%% Prescribed Variables %%
x=-45:.01:45;
dtor = pi/180;
rtod = 180/pi;


%Displacement Vectors:
delta(2) = 0.5+ 4i; %Vector from P1 to P2
delta(3) = 10 + 3i; %Vector from P1 to P3

%Identify PPs
pps(1) = 0+0*1i;
pps(2) = real(delta(2))+imag(delta(2))*1i;
pps(3) = real(delta(3))+imag(delta(3))*1i;


%Set Beta2 and Beta3 (These are only initial values if using loops)
beta2in= 15*dtor; 
beta3in= (20)*dtor;

%Or, instead of using PPs, set Alphas:
alpha(1)=0;
alpha(2)= -45*dtor; %Moving plane rotation from P1 to P2
alpha(3)= -90*dtor; %Moving plane rotatioquivern from P1 to P3

%% Plot Settings
solid = 0; % Circles are either solid or individual points. 1 for solid. 
Von = 0; %Vectors on for Beta3 loops, vectors off (0) forr_ Beta2 loops
firstonly = 0; %1 for only showing vectors on first loop
plotrad = 1; %1 to generate plot of the circle radius for different Beta 2
beta2loops = 180;%3
beta3loops = 180;%2515/2; 

%% Define Pole Locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displacement vector calculation if using PPs
%delta(1) = pps(2) - pps(1); % PP displacement from PP1 to PP2
%delta(2) = pps(3) - pps(1); % PP displacement from PP1 to PP3
%(y-ymid)=m*(x-xmid)

%Start Calculations: 

%Find Pole Locations: 
P12=delta(2)/(1-exp(1i*alpha(2)));
P13=delta(3)/(1-exp(1i*alpha(3)));
P23=((delta(3)*exp(1i*alpha(2)))-(delta(2)*exp(1i*alpha(3))))/(exp(1i*alpha(2))-exp(1i*alpha(3)));
P23pB = -(delta(3)-delta(2))/(exp(1i*alpha(3))-exp(1i*alpha(2)));

%Find Image Pole P23'
yintercept=((real(P13)*imag(P12)-real(P12)*imag(P13))/(real(P13)-real(P12)));
P12P13Slope=((imag(P13)-imag(P12))/(real(P13)-real(P12)));
P12P13Line= x*P12P13Slope+yintercept;
d = (real(P23) + (imag(P23) - yintercept)*P12P13Slope)/(1 + P12P13Slope^2);
P23p=(2*d - real(P23))+1i*(2*d*P12P13Slope - imag(P23) + 2*yintercept);


%% Draw Figure
text = sprintf('Solution Dyads');
figure('Position', [10 50 800 600], 'Name', text,'NumberTitle', 'off'); 
%figure()
hold on;
grid on;
xlim([-10,10]);
ylim([-10,10]);
V=axis;
scale=(V(2)-V(1))/2;

quiver(0,0,1,0,scale,'ok');
quiver(real(delta(2)),imag(delta(2)),real(exp(alpha(2)*1i)),imag(exp(alpha(2)*1i)),scale,'ok'); 
quiver(real(delta(3)),imag(delta(3)),real(exp(alpha(3)*1i)),imag(exp(alpha(3)*1i)),scale,'ok');

plot(real(P23p),imag(P23p),'.','MarkerSize',20,'Color','m')
plot(real(P23),imag(P23),'.','MarkerSize',20,'Color','g')
plot(real(P13),imag(P13),'.','MarkerSize',20,'Color','k')
plot(real(P12),imag(P12),'.','MarkerSize',20,'Color','y')

dataM = zeros(beta3loops,2);
dataK = zeros(beta3loops,2);
r_array = zeros(beta2loops,1);
distance_array = zeros(beta2loops+1,beta3loops/10+1);


%Use distance between poles, interior angle, and trig to solve for radius
%rather than circle system. 

%% Plot
for m=1:beta2loops%1:1:360
    beta2 = (2*m-180)*dtor+beta2in;
    %beta2 = beta2in + 2*pi*(m-1)/beta2loops;
    %beta2 = (m-1)*dtor+beta2in;
    
    for n=1:beta3loops%1:1:360
            beta3 = (2*n-180.1)*dtor; %
            %beta3 = beta3in + 2*pi*(n-1)/beta3loops;%+ .0025*(n-1)*2;
            %beta3 = beta3in + 1*dtor*(n-1);
            %beta3 = (0.5*n-1)*dtor+beta3in;
            
            A = [exp(1i*beta2)-1, exp(1i*alpha(2))-1; exp(1i*beta3)-1, exp(1i*alpha(3))-1];
            B = [delta(2); delta(3)];
            X=A\B;
            W(m,n)=X(1);
            Z(m,n)=X(2);
            
                if(n == 1)
                    distance_array(1,1,m) = beta2*rtod;
                    pminus1 = pps(1)-X(2)-X(1); %On the first pass, set previous as 1st point
                else
                    pcurrent = pps(1)-X(2)-X(1);
                    distance_array(n,2,m) = beta3*rtod;
                    distance_array(n,3,m) = sqrt(real(pcurrent-pminus1)^2+imag(pcurrent-pminus1)^2);
                    pminus1 = pcurrent; %Replace the logged point with the current point
                    
                end

            if(solid == 1)
                dataM(n,:) = [-real(Z(m,n)+W(m,n)),-imag(Z(m,n)+W(m,n))];
                dataK(n,:) = [-real(Z(m,n)),-imag(Z(m,n))];

                if(n == 3)
                  circdatM = findcircle(dataM(1,1),dataM(1,2),dataM(2,1),dataM(2,2),dataM(3,1),dataM(3,2));
                  circdatK = findcircle(dataK(1,1),dataK(1,2),dataK(2,1),dataK(2,2),dataK(3,1),dataK(3,2));
                  %drawcircle(circdatM,1);
                  %drawcircle(circdatK,2);
                end
            
           
            else
                plot(-real(Z(m,n)), -imag(Z(m,n)),'.','Markersize',7, 'Color', 'r'); %Moving
                %plot(-real(Z(m,n))-real(W(m,n)),-imag(Z(m,n))-imag(W(m,n)), '.','Markersize',7,'Color', 'b'); %Ground
            end
            
            if Von==1
               plot([0, -real(Z(m,n))],[0,-imag(Z(m,n))],'r'); %Moving Pivot Circle and coupler are red
               plot([-real(Z(m,n)), -real(Z(m,n))-real(W(m,n))], [-imag(Z(m,n)), -imag(Z(m,n))-imag(W(m,n))],'b');
               %Ground Pivot Circle and links 2/4 are Blue
               if (firstonly == 1)
                   Von = 0;
               end
            end
            %pause(0.025)
            %fprintf('beta 3 is %f\n',beta3*180/pi)
    
    end
    
    if (solid == 1)
    r_array(m,1) = beta2*rtod;
    r_array(m,2) = abs(sqrt(real(P13-P23)^2+imag(P13-P23)^2)/(2*sin(beta2)));
    end
    
    
    %% Plot second direction
    beta3 = beta3in + 2*pi*(m-1)/beta3loops;
    for p=1:beta3loops
        %beta2 = beta2in +2*pi*(p-1)/beta3loops;%+ .0025*(n-1)*2;
        beta2 = (2*p-180)*dtor+beta2in;
        %beta3 = beta3in + 1*dtor*(n-1);
        A = [exp(1i*beta2)-1, exp(1i*alpha(2))-1; exp(1i*beta3)-1, exp(1i*alpha(3))-1];
        B = [delta(2); delta(3)];
%         X=A\B;
        W(m,p)=X(1);
        Z(m,p)=X(2);

        %dataM(p,:) = [-real(Z(m,p)+W(m,p)),-imag(Z(m,p)+W(m,p))];
        dataK(p,:) = [-real(Z(m,p)),-imag(Z(m,p))];
        plot(-real(Z(m,n)), -imag(Z(m,n)),'.','Markersize',7, 'Color', 'r'); %Moving
        %plot(-real(Z(m,n))-real(W(m,n)),-imag(Z(m,n))-imag(W(m,n)), '.','Markersize',7,'Color', 'b'); %Ground
        if(p == 3)
          circdatM = findcircle(dataM(1,1),dataM(1,2),dataM(2,1),dataM(2,2),dataM(3,1),dataM(3,2));
          circdatK = findcircle(dataK(1,1),dataK(1,2),dataK(2,1),dataK(2,2),dataK(3,1),dataK(3,2));
          %drawcircle(circdatM,1);
          %drawcircle(circdatK,2);
        end
    end
end

plot(real(P23p),imag(P23p),'.','MarkerSize',20,'Color','m')
plot(real(P23),imag(P23),'.','MarkerSize',20,'Color','g')
plot(real(P13),imag(P13),'.','MarkerSize',20,'Color','k')
plot(real(P12),imag(P12),'.','MarkerSize',20,'Color','y')

quiver(0,0,1,0,scale,'ok');
quiver(real(delta(2)),imag(delta(2)),real(exp(alpha(2)*1i)),imag(exp(alpha(2)*1i)),scale,'ok'); 
quiver(real(delta(3)),imag(delta(3)),real(exp(alpha(3)*1i)),imag(exp(alpha(3)*1i)),scale,'ok');

quiver(0,0,1,0,scale,'.k');
quiver(real(delta(2)),imag(delta(2)),real(exp(alpha(2)*1i)),imag(exp(alpha(2)*1i)),scale,'.k'); 
quiver(real(delta(3)),imag(delta(3)),real(exp(alpha(3)*1i)),imag(exp(alpha(3)*1i)),scale,'.k');

legend('PP1','PP2','PP3','Image Pole P23','Pole P23', 'Pole P13','Pole P12');
xlim([-15,15]);
ylim([-15,15]);
axis equal;
hold off;


%% Make plots
if(plotrad==1)
    % figure(2)
    % hold on
    % grid on
    % plot(r_array(1:beta2loops,1), r_array(1:beta2loops,2),'Linewidth',2) %Plot showing radius with respect to beta 2
    % xlabel('Beta 2 (Degrees)')
    % ylabel('Radius')
    % title('Circle Radius Relative to Free Choice Beta 2 Value')
    % hold off

    figure(3)
    hold on
    grid on
    colors = jet(9);
    for i = 1:1:(length(r_array)-1)
        if (mod(i,40) == 0)
            gbeta2 = distance_array(1,1,i);
            rnorm = sqrt(real(P13-P23)^2+imag(P13-P23)^2)/(2*sin(gbeta2*dtor));
            a = i/40;
            %sprintf('B2 = %d',gbeta2)
            legstr{a} = sprintf('B2 = %.2f',gbeta2);

            plot(smooth(distance_array(1:beta2loops,2,i)),smooth(distance_array(1:beta2loops,3,i))/abs(rnorm),'Linewidth',2, 'color',colors(i/40,:)) %Plot showing distance between points for all beta 3 at 10 different beta 2s
            %plot(smooth(distance_array(1:359,2,i)),smooth(distance_array(1:359,3,i)),'Linewidth',2,'color',colors(i/40,:)) %Plot showing distance between points for all beta 3 at 10 different beta 2s

        end
    end
    legend(legstr)
    xlabel('Beta 3 (Degrees)')
    ylabel('Distance (normalized by radius)')
    title('Distance between points given Beta 2 and Beta 3')
    hold off
end


%% Functions

function drawcircle(circdat,circtype)
    r = circdat(1);
    center = [real(circdat(2)),imag(circdat(2))];
    pos = [center-r 2*r 2*r];
    if(circtype == 1)
        rectangle('Position',pos,'Curvature',1,'EdgeColor','b','LineWidth',.1)
        %plot(center(1),center(2),'.','MarkerSize',20,'Color','k')
    elseif(circtype == 2)
        rectangle('Position',pos,'Curvature',1,'EdgeColor','r','LineWidth',.1)
        %plot(center(1),center(2),'.','MarkerSize',10,'Color','b')
    end
        %plot(center(1),center(2),'.','MarkerSize',10,'Color','r')
end





