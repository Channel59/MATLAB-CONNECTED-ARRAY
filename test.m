%% Frequency range for analysis
clear 
f0_c   = 17e9; % Design frequency of connected array;
f0_adl = f0_c/2; % Design frequency of ADL; 
f0_eq_permittivity_calc = f0_adl * 1.1; % Permittivity evaluation frequency during synthesis;

f_range = linspace(0.1*f0_c, f0_c);

%% Generate ADL;
Z0 = 80;
ZL = 377;
N = 4;
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
astack = ADLStack(f0_adl, p, er_list);
astack.InitNumMetalLayers = [2,2,2,2];
er_host = 1;
%Synthesize each section
astack.synthesizeIndependently(w_min, w_max, f0_eq_permittivity_calc, AddManufacturingStratification, er_host)
r = astack.getTotalADL();

Zin_ideal = transformer.inputImpedance(r.Frequencies);
Gamma_ideal = (Zin_ideal - Z0) ./ (Zin_ideal + Z0);


%% Genetic algorithm
% 
% numvariables = r.NumMetalLayers * 2; %shift and width;
% 
% x0 = zeros(1,numvariables);
% % Extract initial guess from adl object.
% for i = 1:numvariables / 2
% x0(i) = r.Layers{r.MetalLayerIndices(i)}.w;
% x0(numvariables/2 + i) = r.Layers{r.MetalLayerIndices(i)}.s;
% end
% 
% 
% %constraints: min and max width
% A1 = eye(numvariables);
% A2 = -eye(numvariables);
% A = [A1 ; A2];
% 
% 
% 
% B1a = ones(numvariables/2,1)*w_max;
% B1b = ones(numvariables/2,1)*0.5;
% B2a = -ones(numvariables/2,1)*w_min;
% B2b = ones(numvariables/2,1)*0.5;
% b = [B1a; B1b; B2a; B2b];
% 
% %
% rng default % For reproducibility
% 
% % Objective function
% fun = @(x)simple_objective_10dB_binary(x,r,Gamma_ideal,Z0);
% 
% 
% % Options  
% opts = optimoptions('ga','InitialPopulationMatrix',x0, 'CrossoverFcn', @crossoverintermediate, 'PlotFcn',{'gaplotbestf'} );
% 
% [x,Fval,exitFlag,Output] = ga(fun,numvariables,A,b,[],[],[],[],[],[],opts);
% %%
% for i = 1:numvariables / 2
% r.Layers{r.MetalLayerIndices(i)}.w = x(i);
% r.Layers{r.MetalLayerIndices(i)}.s = x(numvariables/2 + i);
% end
% r.solve()

%%
[Zin_TE, ~] = r.getInputImpedance();
r.Frequencies = f_range;
Gamma_ADL = (Zin_TE - Z0) ./ (Zin_TE + Z0);


%% Connected array

 
l0_c = 3e8 / f0_c;

w = 0.1 * l0_c;
delta_d = 0.1*l0_c;
dx = 0.3*l0_c;
dy = 0.3*l0_c;
h = l0_c /4;
% h = inf;
er = 1;
dielectric_thickness = 0;

ca = ConnectedArray(dx, dy, h, w, delta_d, er, dielectric_thickness);
ca.CavityWalls = true;
ca.ADL_up= r;
c = 0.6*1e-12;
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