function stephenson2_diagram()
% STEPHENSON II six-bar mechanism diagram with adjustable lengths.
% Edit the "params" block below to change geometry. The script solves two
% loop-closure constraints numerically and plots a single pose.

  params = struct();
  % Ground pivots (O2 at origin, O4 along +X)
  params.groundLength = 0.220;        % meters

  % Binary links to the ternary and coupler
  params.inputCrankLength = 0.040;    % O2->P2
  params.couplerLength   = 0.120;     % P2->P3
  params.rockerLength    = 0.110;     % P4->O4
  params.linkO2toP6Len   = 0.090;     % O2->P6

  % Ternary link geometry (three joints: P3, P4, P6)
  params.P3toP4Len       = 0.060;     % P3->P4
  params.P3toP6Len       = 0.080;     % P3->P6
  params.P4P6IncludedDeg = 60.0;      % included angle at P3 between (P3->P4) and (P3->P6)

  % Input crank angle (radians) measured from +X about O2
  params.inputAngleRad   = deg2rad(30);  % adjust as desired

  % Plot controls
  params.showDimensions  = false;
  params.lineWidth       = 2.0;
  params.jointMarkerSize = 60;

  plotStephenson2(params);
end

function plotStephenson2(params)
  O2 = [0.0, 0.0];
  O4 = [params.groundLength, 0.0];

  % Precompute ternary local directions
  beta = deg2rad(params.P4P6IncludedDeg);
  % Place local unit vectors from P3 to P4 and P6 with fixed included angle
  % Reference direction for P3->P4 is angle phi; P3->P6 is at (phi + beta)

  % Known point from input crank
  P2 = O2 + params.inputCrankLength * [cos(params.inputAngleRad), sin(params.inputAngleRad)];

  % Unknowns: theta3 (coupler orientation), phi (ternary orientation)
  % Use fminsearch to minimize squared residuals of two circle constraints
  objective = @(x) closureResidualsSquared(x, O2, O4, P2, params, beta);

  % Initial guess: aim coupler roughly toward O4; ternary roughly horizontal
  theta3_guess = atan2(O4(2) - P2(2), O4(1) - P2(1));
  phi_guess    = 0.0;

  x0 = [theta3_guess, phi_guess];
  options = optimset('Display', 'off');
  x_sol = fminsearch(objective, x0, options);
  theta3 = x_sol(1);
  phi    = x_sol(2);

  % Compute key points
  P3 = P2 + params.couplerLength * [cos(theta3), sin(theta3)];
  e34 = [cos(phi), sin(phi)];
  e36 = [cos(phi + beta), sin(phi + beta)];
  P4 = P3 + params.P3toP4Len * e34;
  P6 = P3 + params.P3toP6Len * e36;

  % Build line segments for plotting (all links black). Ternary drawn as filled triangle below.
  segments = {
    O2, O4;   % ground
    O2, P2;   % input crank
    P2, P3;   % coupler
    P4, O4;   % rocker to right ground
    P6, O2;   % link to left ground
  };

  figure('Color', 'w'); hold on; axis equal;
  for i = 1:size(segments, 1)
    A = segments{i,1};
    B = segments{i,2};
    plot([A(1) B(1)], [A(2) B(2)], '-', 'Color', 'k', 'LineWidth', params.lineWidth);
  end

  % Draw ternary link as filled triangle with pin joints at corners
  V = [P3; P4; P6];
  F = [1 2 3];
  patch('Faces', F, 'Vertices', V, 'FaceColor', 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'k', 'LineWidth', params.lineWidth);

  % Plot joints
  jointPts = [O2; O4; P2; P3; P4; P6];
  scatter(jointPts(:,1), jointPts(:,2), params.jointMarkerSize, 'k', 'filled');


  grid on; box on;
  xlabel('X [m]'); ylabel('Y [m]');
  title('Stephenson II Mechanism (diagram)');

  if params.showDimensions
    annotateDimensions(O2, O4, P2, P3, P4, P6);
  end
end

function err2 = closureResidualsSquared(x, O2, O4, P2, params, beta)
  theta3 = x(1);
  phi    = x(2);

  P3 = P2 + params.couplerLength * [cos(theta3), sin(theta3)];
  e34 = [cos(phi),        sin(phi)       ];
  e36 = [cos(phi + beta), sin(phi + beta)];
  P4 = P3 + params.P3toP4Len * e34;
  P6 = P3 + params.P3toP6Len * e36;

  r1 = (norm(P4 - O4) - params.rockerLength);
  r2 = (norm(P6 - O2) - params.linkO2toP6Len);
  err2 = r1*r1 + r2*r2;
end

function annotateDimensions(O2, O4, P2, P3, P4, P6)
  % Simple dimension ticks for visual reference only.
  dimColor = [0.6 0.6 0.6];
  lw = 1.0;
  plot([O2(1) O4(1)], [O2(2)-0.02 O4(2)-0.02], '-', 'Color', dimColor, 'LineWidth', lw);
  plot([O2(1) O2(1)], [O2(2) O2(2)-0.025], '-', 'Color', dimColor, 'LineWidth', lw);
  plot([O4(1) O4(1)], [O4(2) O4(2)-0.025], '-', 'Color', dimColor, 'LineWidth', lw);
end

% To run: in MATLAB, navigate to this folder and call:
% >> stephenson2_diagram
% Then edit the lengths/angles at the top of this file and re-run.


