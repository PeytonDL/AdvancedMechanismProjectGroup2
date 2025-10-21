%3PP Dyad Motion Synthesis

%% Setup
format short
clear; close;
dtor = (pi / 180); %Degree-to-radian conversion factor
rtod = (180 / pi); %Radian-to-degree conversion factor

%Prescribe Precision Positions:
pps(1) = 0 + 0.0*1i;   % PP1 Initialized at (0,0). 
pps(2) =  4+ 3i;  % PP 2 coords 
pps(3) =  7 + 2i;   % PP 3 coords 
%pps = [0+0i, 4+3i, 7+2i];

% Displacement vector calculation if using Precision Positions
delta(2) = pps(2) - pps(1); % PP displacement from PP1 to PP2
delta(3) = pps(3) - pps(1); % PP displacement from PP1 to PP3

%The delta vectors are what gets used in the standard form equations. Most 
%of the time, from a math perspective, we only care about the displacement
%from one position to the next, not the absolute position of each point.

alpha(1) = 0;
alpha(2) = 30*dtor;
alpha(3) = 65*dtor;% Alpha values are shared by both dyads

%% Solve for Driver Dyad

betaA = [0,-25,-50]*dtor;

betaB = [0,45, 70]*dtor;

output = zeros(1,9);

for j=1:1:1
    

    A=[exp(1i*betaA(2))-1 exp(1i*alpha(2))-1;
       exp(1i*betaA(3))-1 exp(1i*alpha(3))-1];

    deltavect=[delta(2);
               delta(3)];

    Sol_Dyad1=A\deltavect;
    W1D1=Sol_Dyad1(1);
    Z1D1=Sol_Dyad1(2);
    OAD1=pps(1)-Z1D1-W1D1;

    %% Dyad 2 (Output)


    A=[exp(1i*betaB(2))-1 exp(1i*alpha(2))-1;
       exp(1i*betaB(3))-1 exp(1i*alpha(3))-1];

    deltavect=[delta(2);
               delta(3)];

    Sol_Dyad2=A\deltavect;
    W1D2=Sol_Dyad2(1);
    Z1D2=Sol_Dyad2(2);
    OAD2=pps(1)-Z1D2-W1D2;


    %fprintf('D1W: %f+%fi \n',real(W1D1),imag(W1D1))
    %fprintf('D1Z: %f+%fi \n',real(Z1D1),imag(Z1D1))

    %fprintf('D2W: %f+%fi \n',real(W1D2),imag(W1D2))
    %fprintf('D2Z: %f+%fi \n',real(Z1D2),imag(Z1D2))

    %% Calculate transmission angle for present mechanisms
    L1 = norm(W1D1 + Z1D1 - Z1D2 - W1D2);
    L2 = norm(W1D1);
    L3 = norm(Z1D1 - Z1D2);
    L4 = norm(W1D2);
    theta2 = angle(W1D1)-angle(W1D1 + Z1D1 - Z1D2 - W1D2);
    theta4 = angle(W1D2)-angle(W1D1 + Z1D1 - Z1D2 - W1D2);
    theta3 = angle(Z1D1 - Z1D2);
    SideC = norm(OAD1+W1D1-OAD2);
    transmission_ang = acos((L4^2+L3^2-SideC^2)/(2*L4*L3));
    %theta3 = atan2((L1*sin(theta2)-L4*sin(theta4)),(L1*cos(theta2) - L4*cos(theta4)));
    %transmission_ang = atan2(sqrt(1-(((L2^2 + L3^2 - L1^2 - L4^2) - 2*L2*L3*cos(theta2 - theta3))/(2*L2*L3)))^2, ((L2^2 + L3^2 - L1^2 - L4^2) - 2*L2*L3*cos(theta2 - theta3))/(2*L2*L3)); 
    %transmission_ang = abs(theta3-theta4);
    fprintf('Trans Ang: %f \n', transmission_ang*rtod)
    if(transmission_ang*rtod > 90)
        if(transmission_ang*rtod > 180)
            transmission_ang = transmission_ang-pi;
        else
            transmission_ang= pi-transmission_ang;
        end
    end
    %% Store data in output
    output(j,1) = betaA(2);
    output(j,2) = betaA(3);
    output(j,3) = betaB(2);
    output(j,4) = betaB(3);
    output(j,5) = transmission_ang*rtod;
    output(j,6) = W1D1;
    output(j,7) = Z1D1;
    output(j,8) = W1D2;
    output(j,9) = Z1D2;
    %TEST = angle(0+1i)
    
    x = [real(pps(1)), real(pps(1)-Z1D1),real(pps(1)-Z1D2)];
    y = [imag(pps(1)), imag(pps(1)-Z1D1), imag(pps(1)-Z1D2)];
    %z = [real(pps(1)-Z1D2), real(pps(1)-Z1D2), 0];
    c = [.5 .5 .5];
    
    %% Create Figure:
    if(j==1)
        fprintf('Plotting')
        figure(); 
        hold on;
        title('Solution for Dyad Set ',j);
        V=axis;
        scale=(V(2)-V(1))*3;

        for p=1:1:3
            quiver(real(OAD1),imag(OAD1), real(W1D1*exp(1i*betaA(p))),imag(W1D1*exp(1i*betaA(p))),1,'r','linestyle','--','linewidth',2);
            %quiver(real(OAD1+W1D1), imag(OAD1+W1D1),    real(Z1D1),   imag(Z1D1),1,'r');

            quiver(real(OAD2),imag(OAD2),real(W1D2*exp(1i*betaB(p))),imag(W1D2*exp(1i*betaB(p))),1,'b','linestyle',':','linewidth',2);
            %quiver(real(OAD2+W1D2), imag(OAD2+W1D2),    real(Z1D2),   imag(Z1D2),1,'b');

            %plot([real(OAD2+W1D2),real(OAD1+W1D1)],[imag(OAD2+W1D2),imag(OAD1+W1D1)],'k');

            x = [real(pps(p)), real(pps(p)-Z1D1*exp(1i*alpha(p))),real(pps(p)-Z1D2*exp(1i*alpha(p)))];
            y = [imag(pps(p)), imag(pps(p)-Z1D1*exp(1i*alpha(p))), imag(pps(p)-Z1D2*exp(1i*alpha(p)))];
            c = [.5 .5 .5];

            fill(x,y,c,'FaceAlpha',0.5);

            %h.FaceAlpha = 0.5;
            quiver(real(pps(1)),imag(pps(1)),real(exp(alpha(1)*1i)),imag(exp(alpha(1)*1i)),scale,'ok');
            quiver(real(delta(2))+real(pps(1)),imag(delta(2))+imag(pps(1)),real(exp(alpha(2)*1i)),imag(exp(alpha(2)*1i)),scale,'xk'); 
            quiver(real(delta(3))+real(pps(1)),imag(delta(3))+imag(pps(1)),real(exp(alpha(3)*1i)),imag(exp(alpha(3)*1i)),scale,'dk');

            plot(real(pps(1)-Z1D1-W1D1),imag(pps(1)-Z1D1-W1D1),'*r');
            plot(real(pps(1)-Z1D2-W1D2),imag(pps(1)-Z1D2-W1D2),'*b');


            
            
        end
    end
end
legend('Input Link', 'Follower Link', 'Coupler', 'PP1','PP2','PP3','OA','OB')
axis equal
hold off

