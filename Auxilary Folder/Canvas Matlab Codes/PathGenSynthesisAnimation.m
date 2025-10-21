%3PP Dyad Motion Synthesis

%% Setup
format short
clear; close;
dtor = (pi / 180); %Degree-to-radian conversion factor
rtod = (180 / pi); %Radian-to-degree conversion factor

%Prescribe Precision Positions:
pps(1,1) = 0 + 0.0*1i;   % PP1 Initialized at (0,0). 
pps(2,1) =  -2+ 3i;  % PP 2 coords 
pps(3,1) =  -1 + 6i;   % PP 3 coords 

% Displacement vector calculation if using Precision Positions
delta(2) = pps(2) - pps(1); % PP displacement from PP1 to PP2
delta(3) = pps(3) - pps(1); % PP displacement from PP1 to PP3

%The delta vectors are what gets used in the standard form equations. Most 
%of the time, from a math perspective, we only care about the displacement
%from one position to the next, not the absolute position of each point.

pps(1,2) = 0; %Alpha angles as second column of pps matrix.
pps(2,2) = -30*dtor;
pps(3,2) = -75*dtor;% Alpha values are shared by both dyads

%% Solve for Driver Dyad

beta2A = 35*dtor;%[-25 60 90 30 30 30]*dtor;
beta3A = 85*dtor;

beta2B = -10*dtor;%[45 20 20 -20 -50 -70]*dtor;
beta3B = -55*dtor;

output = zeros(length(beta2A),5);

for j=1:1:length(beta2A)
    

    A=[exp(1i*beta2A(j))-1 exp(1i*pps(2,2))-1;
       exp(1i*beta3A)-1 exp(1i*pps(3,2))-1];

    deltavect=[delta(2);
               delta(3)];

    Sol_Dyad1=A\deltavect;
    W1D1=Sol_Dyad1(1);
    Z1D1=Sol_Dyad1(2);
    OAD1=pps(1)-Z1D1-W1D1;

    %% Dyad 2 (Output)

    A=[exp(1i*beta2B(j))-1 exp(1i*pps(2,2))-1;
       exp(1i*beta3B)-1 exp(1i*pps(3,2))-1];

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


    output(j,1) = beta2A(j);
    output(j,2) = beta3A;
    output(j,3) = beta2B(j);
    output(j,4) = beta3B;
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
    
end


Link1 = W1D1 + Z1D1 - Z1D2 - W1D2;
Link2 = W1D1;
Link3 = Z1D1 - Z1D2;
Link4 = W1D2;

inZang = angle(Z1D1)-angle(Link3);

%% Analyze Grashof Type
linkagemag = [abs(Link1), abs(Link2), abs(Link3), abs(Link4)];
grashcheck = sort(linkagemag);


if((grashcheck(1)+grashcheck(4))<=(grashcheck(2)+grashcheck(3)))
    disp('This is a Grashof Mechanism')
    grashof = true;
    if(grashcheck(1)==abs(Link2))
        disp('Driver link is shortest')
    end
else
    disp('This is not a Grashof Mechanism')
    grashof = false;
end

%% Animation Parameters and Run Animation
num_frames = 360;
frameduration = 0.001; %Decrease for faster movement
direction = -1; %+1 for CCW Input rotation, -1 for CW Input rotation
theta_step = 2*pi/num_frames; % Step size for input rotation
input_link = 1; % Choose 1 or 4 for the input link

four_bar_animation(Link1, Link2, Link3, Link4, Z1D1, pps, theta_step,num_frames, frameduration, direction);

% MATLAB Code: Animation of a 4-bar linkage mechanism
% Define the lengths of the links as complex numbers
% Input: Link lengths (L1, L2, L3, L4), initial positions, input angle


%% Animation Function %%

function four_bar_animation(L1, L2, L3, L4, Z1, pps, theta_step,num_frames, frameduration, direction)
    % Parameters:
    % L1, L2, L3, L4: lengths of the links as complex numbers
    % theta_start: initial angle of the input link (in radians)
    % theta_step: step size for input link rotation (in radians)
    % input_link: specify which link (1-4) is the input link

    % Initialize
    theta = 0; % Start angle for the input link
    inZang = angle(Z1)-angle(L3);
    
    delta(2) = pps(2) - pps(1); % PP displacement from PP1 to PP2
    delta(3) = pps(3) - pps(1); % PP displacement from PP1 to PP3

    % Set up the animation figure
    figure(1);
    clf;
    axis equal;
    grid on;
    hold on;
    title('4-Bar Linkage Mechanism Animation');
    xlabel('X'); 
    ylabel('Y');

    % Base positions
    A0 = pps(1)-Z1-L2; % Fixed base joint (ground)
    B0 = A0+L1; % Other fixed joint (ground)
    max_extent = max([abs(L1), abs(L2), abs(L3), abs(L4)]) * 2; % Adjust multiplier if necessary
    axis([-max_extent, max_extent, -max_extent/2, max_extent+max_extent/2]); % Set fixed axis limits

    
    % Choose the valid solution for the coupler
    linkageposition = zeros(num_frames,6);
    
    %% Gif Parameters
    Th2Res = 100;
    SecsPerBranch = 3;
    FPS = Th2Res/SecsPerBranch; % Frame Rate
    fileName = sprintf('PathGenAnimation2.gif');
    
    %% Make Plots
    
    for i=1:num_frames
    % Plot the linkage
    theta = theta + direction*theta_step;
    [L2new, L3new, L4new, Znew, transang, thetao] = analytical_displacement(L1, L2, L3, L4, Z1, inZang, A0,theta);
    linkageposition(i,:) = [L2new, L3new, L4new, Znew, thetao, transang];
    scale = 2.5;
    end
    mintransang = min(linkageposition(:,6));
    fprintf('Minimum transmission angle is %f', mintransang*180/pi);
    for j=1:length(linkageposition(:,1))
    
    cla;
    %Plot Prescribed Positions
    quiver(real(pps(1)),imag(pps(1)),1,0,scale,'ok');
    quiver(real(pps(1)+delta(2)),imag(pps(1)+delta(2)),real(exp(pps(2,2)*1i)),imag(exp(pps(2,2)*1i)),scale,'ob'); 
    quiver(real(pps(1)+delta(3)),imag(pps(1)+delta(3)),real(exp(pps(3,2)*1i)),imag(exp(pps(3,2)*1i)),scale,'or');
   
    
    %Plot Links
    plot([real(A0), real(A0+linkageposition(j,1))], [imag(A0), imag(A0+linkageposition(j,1))], 'r', 'LineWidth', 2); % Link 2
    plot([real(A0+linkageposition(j,1)), real(A0+linkageposition(j,1)+linkageposition(j,2))], [imag(A0+linkageposition(j,1)), imag(A0+linkageposition(j,1)+linkageposition(j,2))], 'g', 'LineWidth', 2); % Link 3
    plot([real(B0), real(B0+linkageposition(j,3))], [imag(B0), imag(B0+linkageposition(j,3))], 'b', 'LineWidth', 2); % Link 4
    %Plot Coupler Triangle
    plot([real(A0+linkageposition(j,1)), real(A0+linkageposition(j,1)+linkageposition(j,4))], [imag(A0+linkageposition(j,1)), imag(A0+linkageposition(j,1)+linkageposition(j,4))], 'g', 'LineWidth', 2);% Link 3 Side
    plot([real(A0+linkageposition(j,1)+linkageposition(j,4)),real(B0+linkageposition(j,3))], [imag(A0+linkageposition(j,1)+linkageposition(j,4)),imag(B0+linkageposition(j,3))], 'g', 'LineWidth', 2);% Link 3 Side
    
    % Plot the joints
     plot(real(A0), imag(A0), 'ko', 'MarkerFaceColor', 'k','MarkerSize',5); % A0
     plot(real(A0+linkageposition(j,1)), imag(A0+linkageposition(j,1)), 'ko', 'MarkerFaceColor', 'k','MarkerSize',3); % A
     plot(real(B0), imag(B0), 'ko', 'MarkerFaceColor', 'k','MarkerSize',5); % B0
     plot(real(B0+linkageposition(j,3)), imag(B0+linkageposition(j,3)), 'ko', 'MarkerFaceColor', 'k','MarkerSize',3); % B

    
    % Get frame, convert to image, index colors
    frame = getframe(gcf);
    im = frame2im(frame);
    [indexedImage,Colormap] = rgb2ind(im,256);
    % if first frame, create GIF file
    if j == 1
        imwrite(indexedImage,Colormap,fileName,'gif','LoopCount',Inf,'DelayTime',1/FPS);
    % else, append gif with this frame
    else
        imwrite(indexedImage,Colormap,fileName,'gif','WriteMode', 'append','DelayTime',1/FPS);
    end


    pause(frameduration); % Pause for animation effect
    end
    hold off;
end



%% Analytical Displacement Function

    function [L2new, L3new, L4new, Znew, transang, theta] = analytical_displacement(L1, L2, L3, L4, Z1, ~, ~, theta)
    %B0 = A0 + L1;
    
    %Find r7*
    r7i = L2-L1;
    %Find angle between r7 and link 4
    psi = (angle(L4)-angle(r7i))*(180 / pi);
    if(psi > 180)
        psi = psi - 360;
    elseif(psi < -180)
        psi = psi + 360;
    end
    %Find sign of mu
    if(psi>0)
        mu = 1;
    elseif(psi<0)
        mu = -1;
    else
        disp('Change Point Configuration pos 1')
    end
    
    %(4) Calculate new L2
    L2new = L2 * exp(1i * theta);
    
    %(5) Find r7 new
    r7 = L2new - L1;
    
    %Calculate new psi: 
    psinew = acos((abs(L4)^2 + abs(r7)^2 - abs(L3)^2)/(2*abs(L4)*abs(r7)));
    
    %(7) Calculate new Link 4
    L4new = abs(L4)* exp(1i * (angle(r7) + psinew*mu));
    
    %(8)
    L3new = L4new-r7;
    %magL3new = abs(L3new);
    %Find new Z vector for coupler vertex
    Znew = Z1*exp(1i*(angle(L3new) -angle(L3)));
    
    % Transmission Angle Calculations
    
    transang = acos((abs(L4)^2+abs(L3)^2-abs(r7)^2)/(2*abs(L4)*abs(L3)));
    %fprintf('Trans Ang: %f \n', transang*rtod)
    if(transang*pi/180 > 90)
        if(transang*rtod > 180)
            transang = transang-pi;
        else
            transang= pi-transang;
        end
    
    end
end
    
