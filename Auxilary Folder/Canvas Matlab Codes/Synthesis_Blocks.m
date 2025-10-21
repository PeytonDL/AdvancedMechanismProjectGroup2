%% Synthesis Blocks Program
% Example Created for ME 5243
% 2/26/24
% In multiloop mechanism synthesis, one way to approach the problem is to
% divide the intended topology into a collection of discrete, dependent
% chains. These chains can  be considered as the "building blocks" of complex 
%mechanisms. This program offers examples of several types of the base
%building blocks that designers can use to assemble larger mechanisms.
%% Setup
clear; csc; close;
dtor = (pi / 180); %Degree-to-radian conversion factor
rtod = (180 / pi); %Radian-to-degree conversion factor

%Prescribe Precision Positions:
pps(1) = 0 + 0.0*1i;   % PP1 Initialized at (0,0). 
pps(2) =  1+ 2i;  % PP 2 coords 
pps(3) =  2 + 1i;   % PP 3 coords 
pps(4) = 3 + 2i;   % PP 4 coords 

% Displacement vector calculation if using Precision Positions
delta(2) = pps(2) - pps(1); % PP displacement from PP1 to PP2
delta(3) = pps(3) - pps(1); % PP displacement from PP1 to PP3
delta(4) = pps(4) - pps(1); % PP displacement from PP1 to PP2
%The delta vectors are what gets used in the standard form equations. Most 
%of the time, from a math perspective, we only care about the displacement
%from one position to the next, not the absolute position of each point.

%For motion generation problems, we prescribe a position and an angle
%at each precision position. Here gamma refers to the angle of the third
%link in the triad. 
gamma(1) = 0;
gamma(2) = 5*dtor;
gamma(3) = 50*dtor;
gamma(4) = -25*dtor;

alpha(1) = 0;
alpha(2) = 5*dtor;
alpha(3) = 50*dtor;% Alpha values are shared by both dyads
alpha(4) = 75*dtor;
%The "dtor" operator is used to convert from degrees to radians. Degrees
%tend to be more intuitive for users, but in the majority of the program
%it is a good idea to operate in radians to avoid errors.

%% General Plotting
figure() %Increment figure number to keep previous figures open
hold on
text = sprintf('Solution Chains');
%figure('Position', [10 50 800 600], 'Name', text,'NumberTitle', 'off'); 
grid on;
title('Here is a figure');
h=axis;
scale=(h(2)-h(1))/5;
%Plot precision positions: 
quiver(real(pps(1)),imag(pps(1)),1,0,scale,'ok');
quiver(real(pps(1))+real(delta(2)),imag(pps(1))+imag(delta(2)),real(exp(gamma(2)*1i)),imag(exp(gamma(2)*1i)),scale,'ok'); 
quiver(real(pps(1))+real(delta(3)),imag(pps(1))+imag(delta(3)),real(exp(gamma(3)*1i)),imag(exp(gamma(3)*1i)),scale,'ok');
quiver(real(pps(1))+real(delta(4)),imag(pps(1))+imag(delta(4)),real(exp(gamma(4)*1i)),imag(exp(gamma(4)*1i)),scale,'ok');

axis equal;
%hold off; %Make sure to end each unique figure with "hold off" to stop
%adding to the same figure. 

%% Dyad in Three Prescribed Positions (Motion Generation)
betaD1(2) = -25*dtor; %Initial choice for Beta 2, converted to radians
betaD1(3) = -55*dtor;


alpha(1) = 0;
alpha(2) = 5*dtor;
alpha(3) = 50*dtor;% Alpha values are shared by both dyads

A=[exp(1i*betaD1(2))-1 exp(1i*alpha(2))-1;
   exp(1i*betaD1(3))-1 exp(1i*alpha(3))-1];

deltavect=[delta(2);
           delta(3)];

Sol_Dyad1=A\deltavect;
W1D1=Sol_Dyad1(1);
Z1D1=Sol_Dyad1(2);
OAD1=pps(1)-Z1D1-W1D1;

%% Dyad in Four Prescribed Positions (Motion Generation)(Compatibility Linkage Approach)
A2_2 = 6*dtor;
beta(1)=0*dtor; 
beta(2)= -19*dtor; %Moving plane rotation from P1 to P2
beta(3)= -21*dtor; %Moving plane rotation from P1 to P3
beta(4)= -12*dtor; %Moving plane rotation from P1 to P4
% Set Free choice alpha2


%Set up Compatibility Equations: 
%Define the delta (uppercase/triangle) values
d2=det([exp(1i*beta(3))-1 delta(3); 
    exp(1i*beta(4))-1 delta(4)]);
d3=(-1*det([exp(1i*beta(2))-1 delta(2); 
    exp(1i*beta(4))-1 delta(4)]));
d4=(det([exp(1i*beta(2))-1 delta(2); 
    exp(1i*beta(3))-1 delta(3)]));
d1=(-d2-d3-d4);

%Solve for B3 and B4
d=d1+d2*exp(1i*A2_2);
x=(abs(d4)^2-abs(d3)^2-abs(d)^2)/(2*abs(d3)*abs(d));
A3_2=(angle(d)+(acos(x))-angle(d3)); %Issue with Acos(x) domain here
x2=(abs(d3)^2-abs(d4)^2-abs(d)^2)/(2*abs(d4)*abs(d));
A4_2= (angle(d) - acos(x2) - angle(d4));
if(x>1 || x<-1)
    C='Bad Solution Dyad One'
    %Bug fixing step required here to root out solutions with abs(x)>1
end
A=([exp(1i*A2_2)-1 exp(1i*beta(2))-1;
   exp(1i*A3_2)-1 exp(1i*beta(3))-1]);
deltavect=[delta(2);
           delta(3)];

Sol_Dyad1=A\deltavect;
Z1D2=Sol_Dyad1(1);
W1D2=Sol_Dyad1(2);
OAD2=pps(1)-Z1D2-W1D2;
alpha2=beta(2);
alpha3=beta(3);%B3_2;
alpha4=beta(4);%B4_2;
%% Triad in Three Prescribed Positions (Motion Gen + Ground Pivot Specification)

R(1)= X + Yi; %R1 extends from O to pps1
R(2)=delta(2)+R(1);
R(3)=delta(3)+R(1);


A = [1 1 1;
    exp(1i*beta(2)), exp(1i*alpha(2)), exp(1i*gamma(2));
    exp(1i*beta(3)), exp(1i*alpha(3)), exp(1i*gamma(3))]; %Added alpha4-gamma4, remove for linear solution

Rvect=[R(1);
       R(2);
       R(3)]; %Added R(4), take away for linear solution. 

X=A\Rvect;
W1T=X(1);
Z1T=X(2);
V1T=X(3);
%% Triad in Four Prescribed Positions (Motion Gen with Prescribed Timing)
A= [exp(1i*beta(2))-1, exp(1i*alpha(2))-1, exp(1i*gamma(2))-1; 
        exp(1i*beta(3))-1, exp(1i*alpha(3))-1, exp(1i*gamma(3))-1;
        exp(1i*beta(4))-1, exp(1i*alpha(4))-1, exp(1i*gamma(4))-1];
            
    B=[delta(2);
       delta(3);
       delta(4)];

    X=A\B;
    W =X(1);
    V =X(2);
    Z =X(3);
%% Triad in Four Prescribed Positions with Ground Pivot Specification (Motion Gen)
R(1)=X + Yi; 
R(2)=delta(2)+R(1);
R(3)=delta(3)+R(1);
R(4)=delta(4)+R(1);


d2=det([exp(1i*alpha2) exp(1i*gamma(2)) R(2); 
    exp(1i*alpha3) exp(1i*gamma(3)) R(3);
    exp(1i*alpha4) exp(1i*gamma(4)) R(4)]);

d3=-1*det([1 1 R(1); 
    exp(1i*alpha3) exp(1i*gamma(3)) R(3);
    exp(1i*alpha4) exp(1i*gamma(4)) R(4)]);

d4=det([1 1 R(1); 
    exp(1i*alpha2) exp(1i*gamma(2)) R(2);
    exp(1i*alpha4) exp(1i*gamma(4)) R(4)]);

d5=-1*det([1 1 R(1); 
    exp(1i*alpha2) exp(1i*gamma(2)) R(2);
    exp(1i*alpha3) exp(1i*gamma(3)) R(3)]);

d=d2+d3*exp(1i*beta(2));%+d1;

x=((abs(d5)^2-abs(d4)^2-abs(d)^2)/(2*abs(d4)*abs(d)));
B3T=(angle(d)+(acos(x))-angle(d4));
x2=(abs(d4)^2-abs(d5)^2-abs(d)^2)/(2*abs(d5)*abs(d));
if(x>1 || x<-1)
    x
    C='Bad Solution Triad'
    %Bug fixing step here
end
B4T=(angle(d)-acos(x2)-angle(d5));

A = [1 1 1;
    exp(1i*alpha2), exp(1i*beta(2)), exp(1i*gamma(2));
    exp(1i*alpha3), exp(1i*B3T), exp(1i*gamma(3))]; 

Rvect=[R(1);
       R(2);
       R(3)]; 

X=A\Rvect;
W1T=X(1);
V1T=X(2);
Z1T=X(3);
beta(3)=B3T;
beta(4)=B4T;
Z2=Z1T*exp(1i*gamma(2));
Z3=Z1T*exp(1i*gamma(3));
Z4=Z1T*exp(1i*gamma(4));
OAT=(real(pps(1))-real(Z1T)-real(V1T)-real(W1T))+1i*(imag(pps(1))-imag(Z1T)-imag(V1T)-imag(W1T));
%Triad (one) Data
alphaT = [0 alpha2 alpha3 alpha4];
betaT = [0 beta(2) beta(3) beta(4)];
gammaT = [0 gamma(2) gamma(3) gamma(4)];
%% End Code
%CSV Output if desired: 
output = [1, 2, 3, 4;
          5, 6, 7, 8];
%csvwrite('This_Is_An_Output_File.csv',output)

%The above line outputs a csv file to your working directory with title
%"This_Is_An_Output_File.csv" If you change some variables and run the program again without
%changing the title "output" here, the program will automatically overwrite
%the previous file. Useful if the data wasn't any good, but frustrating if
%you wanted to hang on to a particular result.

hold off; %Stop adding to the current figure. 