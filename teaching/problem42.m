% ==========================================================================================
% [Solid Mechanics] Problem 42: Analysis of torsion of a prismatic shaft with triangular
%                               section using Saint-Venant Theory
%
% António Manuel Couto Carneiro < amcc@fe.up.pt >
% Department of Mechanical Engineering @ FEUP, 2020
% ==========================================================================================
clearvars
%% Initialise variables and define the Saint-Venant function
syms x y h Phi K G theta ... % Problem variables
     C1 C2 C3 ...            % Implicit definition of the boundaries
     y0 y1 x0 x1 ...         % Auxiliary variables for defining the boundaries
     xC1 xC2 xC3 ...         % Boundary equations expressed as a function of y only
     valC1 valC2 valC3 ...   % Value of the potential function at the boundaries
     yC1 yC2 ...             % Boundary equations expressed as a function of x only
     tauxz tauyz             % Shear stress componenents

Phi = K * (x - sqrt(3) * y - 2 * h / 3) * (x + sqrt(3) * y - 2 * h / 3) * (x + h / 3);

fprintf("Saint-Venant function")
fprintf("\n=====================\n")
pretty(Phi)

%% Expand the function
Phi = expand(Phi);
fprintf("\nSaint-Venant function (symbolic expansion)")
fprintf("\n==========================================\n")
pretty(Phi)

%% Define the boundaries
% Bottom side
x0 = -h / 3;
x1 = 2 * h / 3;
y0 = -h / sqrt(3);
y1 = 0;
C1 = y - y0 - (y1 - y0) / (x1 - x0) * (x - x0);

% Top side
x0 = -h / 3;
x1 = 2 * h / 3;
y0 = h / sqrt(3);
y1 = 0;
C2 = y - y0 - (y1 - y0) / (x1 - x0) * (x - x0);

% Left side
C3 = x + h / 3;

fprintf("\nBottom boundary")
fprintf("\n================\n")
pretty(collect(C1, x))

fprintf("\nTop boundary")
fprintf("\n================\n")
pretty(collect(C2, x))

fprintf("\nLeft boundary")
fprintf("\n================\n")
pretty(collect(C3, x))

%% Check if the value of the potential function is zero at the boundaries
% Bottom side
xC1 = solve(C1 == 0, x);
valC1 = subs(Phi, x, xC1);
valC1 = simplify(valC1);

% Top side
xC2 = solve(C2 == 0, x);
valC2 = subs(Phi, x, xC2);
valC2 = simplify(valC2);

% Left side
xC3 = solve(C3 == 0, x);
valC3 = subs(Phi, x, xC3);
valC3 = simplify(valC3);

fprintf("\nValue at the bottom boundary")
fprintf("\n============================\n")
pretty(valC1)

fprintf("\nValue at the top boundary")
fprintf("\n=========================\n")
pretty(valC2)


fprintf("\nValue at the side boundary")
fprintf("\n==========================\n")
pretty(valC3)


%% Perform partial differentiation
% First and second derivatives in order to x
dPhidx = diff(Phi, x);
d2Phidx2 = diff(dPhidx, x);

% First and second derivaties in order to y
dPhidy = diff(Phi, y);
d2Phidy2 = diff(dPhidy, y);

fprintf("\nFirst partial derivative of Phi in order to x: dPhi / dx")
fprintf("\n========================================================\n")
pretty(dPhidx)

fprintf("\nFirst partial derivative of Phi in order to y: dPhi / dy")
fprintf("\n========================================================\n")
pretty(dPhidy)

fprintf("\nSecond partial derivative of Phi in order to x: d2Phi / dx2")
fprintf("\n===========================================================\n")
pretty(d2Phidx2)

fprintf("\nSecond partial derivative of Phi in order to y: d2Phi / dy2")
fprintf("\n===========================================================\n")
pretty(d2Phidy2)

%% Determine the constant K in order to verify the compatibility equation
equation = d2Phidx2 + d2Phidy2 == -2 * G * theta;
Ksol = solve(equation, K);

fprintf("\nConstant K satisfying the compatibility equation")
fprintf("\n================================================\n")
pretty(Ksol)

%% Compute the torsional moment
% Integration limits in the y direction
yC1 = simplify(solve(C1 == 0, y));
yC2 = simplify(solve(C2 == 0, y));

% Integrate and substitute the value of K
Mt = 2 * int(int(Phi, y, [yC1 yC2]), x, [-h / 3,  2 * h / 3]);
Mt = subs(Mt, K, Ksol);

fprintf("\nTorsional Moment")
fprintf("\n================\n")
pretty(Mt)

%% Determine the shear stress components
tauxz = subs(-dPhidy, K, Ksol);
tauyz = subs(dPhidx, K, Ksol);

fprintf("\nXZ Shear stress component")
fprintf("\n=========================\n")
pretty(simplify(tauxz))

fprintf("\nYZ Shear stress component")
fprintf("\n=========================\n")
pretty(simplify(tauyz))
% ==========================================================================================