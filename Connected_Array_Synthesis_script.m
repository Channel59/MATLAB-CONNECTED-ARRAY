%% Frequency range for analysis
clear 
f0_c   = 32e9; % Design frequency of connected array;
f0_adl = 22e9; % Design frequency of ADL; 
f0_eq_permittivity_calc = f0_adl * 1.1; % Permittivity evaluation frequency during synthesis;

f_range = linspace(0.1*f0_c, f0_c,20);
% f_range = 1e9;
%% Generate ADL;
Z0 = 80;
ZL = 377;
N = 2;
Gamma_m = 0.05;
l0_adl = 3e8 / f0_adl;
p = 0.1*l0_adl;
AddManufacturingStratification = false;


w_min = 0.1e-3;
w_max = 0.7*p;

% Transformer
transformer = ChebyshevTransformer(Z0, ZL, N, Gamma_m, f0_adl);



% Find required permittivities
er_list = (120*pi  ./ transformer.impedances()) .^2;

% Create ADL stack
astack = ADLStack(f0_adl, p, er_list, f_range);
astack.InitNumMetalLayers = [2,2,2,2];
er_host = 1;
%Synthesize each section
astack.synthesizeIndependently(w_min, w_max, f0_eq_permittivity_calc, AddManufacturingStratification, er_host)
r = astack.getTotalADL();
r.Frequencies = f_range;
r.solve();
Zin_ideal = transformer.inputImpedance(f_range);
Gamma_ideal = (Zin_ideal - Z0) ./ (Zin_ideal + Z0);


%% Genetic algorithm

numvariables = r.NumMetalLayers * 2; %shift and width;

x0 = zeros(1,numvariables);
% Extract initial guess from adl object.
for i = 1:numvariables / 2
x0(i) = r.Layers{r.MetalLayerIndices(i)}.w;
x0(numvariables/2 + i) = r.Layers{r.MetalLayerIndices(i)}.s;
end


%constraints: min and max width
A1 = eye(numvariables);
A2 = -eye(numvariables);
A = [A1 ; A2];



B1a = ones(numvariables/2,1)*w_max;
B1b = ones(numvariables/2,1)*0.5;
B2a = -ones(numvariables/2,1)*w_min;
B2b = ones(numvariables/2,1)*0.5;
b = [B1a; B1b; B2a; B2b];

%
rng default % For reproducibility

% Objective function
fun = @(x)simple_objective(x,r,Gamma_ideal,Z0);


% Options  
opts = optimoptions('ga','InitialPopulationMatrix',x0, 'CrossoverFcn', @crossoverintermediate, 'PlotFcn',{'gaplotbestf'} );

[x,Fval,exitFlag,Output] = ga(fun,numvariables,A,b,[],[],[],[],[],[],opts);
%
for i = 1:numvariables / 2
r.Layers{r.MetalLayerIndices(i)}.w = x(i);
r.Layers{r.MetalLayerIndices(i)}.s = x(numvariables/2 + i);
end
r.solve()

%%
r.Frequencies = f_range;
[Zin_TE, ~] = r.getInputImpedance(f_range, eps);

Gamma_ADL = (Zin_TE - Z0) ./ (Zin_TE + Z0);


%% Connected array

 
l0_c = 3e8 / f0_c;

w = 0.1 * l0_c;
delta_d = 0.1*l0_c;
dx = p;
dy = p;
h = l0_c /4;
% h = inf;
er = 2.2;
dielectric_thickness = 0;

ca = ConnectedArray_new(dx, dy, h, w, delta_d, er, dielectric_thickness);
ca.CavityWalls = true;
ca.ADL_up= r;
c = 11*1e-12;
C = 1./(1j * 2*pi *f_range*c);
% C = 0;
Z_in = ca.activeInputImpedance(f_range, deg2rad(eps),deg2rad(90)) + C;
% Z_in = ca.activeInputImpedance(4e9, deg2rad(eps),eps) + C;
% 
figure;
subplot(1,2,1)
plot(f_range, real(Z_in))
hold on
plot(f_range, imag(Z_in))
grid on
ylim([-50 200])


subplot(1,2,2)
Gamma = (Z_in - Z0) ./ (Z_in + Z0);
plot(f_range, 20*log10(abs(Gamma)), 'LineWidth', 2)
hold on 
plot(f_range, 20*log10(abs(Gamma_ADL)), 'LineWidth', 2)
grid on
legend("|\Gamma| (Connected array)","|\Gamma| (ADL)")


%% Connected array genetic algorithm
%%  variables considered:
%   C - Capacitance, min value: 0, max value: 1e-10 Farad;
Cmin = 0;
Cmax = 1e-10;
%   w - slot width
wmin = 0.8e-3;
wmax = 0.5*dx;
%   delta - feed width
delta_min = 0.1*dx;
delta_max = 0.6*dx;
%   h - backing reflector height
h_min = l0_c/15;
h_max = l0_c/4;



numvariables = 4;

x0 = zeros(1,numvariables);
% Extract initial guess from adl object.
x0(1) = c;
x0(2) = ca.w;
x0(3) = ca.delta_d;
x0(4) = ca.h;


%constraints: min and max width
A1 = eye(numvariables); % maximum values 
A2 = -eye(numvariables); % minimum values
A = [A1 ; A2];
B1 = [Cmax; wmax; delta_max; h_max]; % Maximum values
B2 = [Cmin; wmin; delta_min; h_min]; % Minimum values
b = [B1; -B2];

%
rng default % For reproducibility
%
% Pre-compute ADL input impedance;
[Z_up_TE, Z_up_TM] = ca.getADLimpedance(f_range, deg2rad(eps) ,deg2rad(eps));

% Objective function


fun = @(x)connected_array_objective(x,ca,Gamma_ideal,Z0, f_range, Z_up_TE, Z_up_TM);


% Options  
opts = optimoptions('ga','InitialPopulationMatrix',x0, 'CrossoverFcn', @crossoverintermediate, 'PlotFcn',{'gaplotbestf'} );

[x,Fval,exitFlag,Output] = ga(fun,numvariables,A,b,[],[],[],[],[],[],opts);
% Save the solution
c = x(1);
ca.w = x(2);
ca.delta_d = x(3);
ca.h = x(4);

fprintf("Design finished -- Summary\n")
fprintf("del:\t%f mm \t%f l0\n", ca.delta_d*1e3, ca.delta_d/l0_c)
fprintf("w:\t\t%f mm \t%f l0\n", ca.w*1e3, ca.w/l0_c)
fprintf("h:\t\t%f mm \t%f l0\n", ca.h*1e3 ,ca.h/l0_c)
fprintf("dx:\t\t%f mm \t%f l0\n", ca.dx*1e3 ,ca.dx/l0_c)
fprintf("C:\t\t%f pF\n", c*1e12)
fprintf("p:\t\t%f mm \t%f l0\n", p*1e3 ,p/l0_c)


% C = 0;
%%
Theta = deg2rad(eps);
Phi = deg2rad(90);

f_range = linspace(0.1*f0_c, f0_c,100);
C = 1./(1j * 2*pi *f_range*c);
Zin_ideal = transformer.inputImpedance(f_range);
Gamma_ideal = (Zin_ideal - Z0) ./ (Zin_ideal + Z0);
r.Frequencies = f_range;
[Zin_TE, ~] = r.getInputImpedance(f_range, Theta);

Gamma_ADL = (Zin_TE - Z0) ./ (Zin_TE + Z0);
Z_in = ca.activeInputImpedance(f_range, Theta,Phi) + C;
% Z_in = ca.activeInputImpedance(4e9, deg2rad(eps),eps) + C;
% 
figure(15);
clf('reset') 
subplot(1,2,1)
plot(f_range, real(Z_in))
hold on
plot(f_range, imag(Z_in))
grid on
ylim([-50 200])


subplot(1,2,2)
Gamma = (Z_in - Z0) ./ (Z_in + Z0);
plot(f_range, 20*log10(abs(Gamma)), 'LineWidth', 2)
hold on 
plot(f_range, 20*log10(abs(Gamma_ADL)), 'LineWidth', 2)
grid on
legend("|\Gamma| (Connected array)","|\Gamma| (ADL)")
title("\theta = " + rad2deg(Theta))


