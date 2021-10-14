classdef ConnectedArray_new < handle
    %CONNECTEDARRAY Class that implements a connected array
    %   Detailed explanation goes here
    
    properties
        f;
        dy;
        dx;
        h; %reflector height
        w % slot/dipole thickness
        er_back;
        dielectric_thickness; % back dielectric thickness;
        delta_d;
        theta;
        phi;
        type; % slot or strip
        ADL_up;
        CavityWalls = false;
    end
    
    methods
        function obj = ConnectedArray_new(dx, dy, h, w, delta_d,er, dielectric_thickness)
            %CONNECTEDARRAY Construct an instance of this class
            %   Detailed explanation goes here
            obj.dx = dx;
            obj.dy = dy;
            obj.h  = h;
            obj.w = w;
            obj.delta_d = delta_d;
            obj.dielectric_thickness = dielectric_thickness;
            obj.er_back = er;
            obj.ADL_up = ADL();
        end
        
        function Z_a = activeInputImpedance(obj, f, theta ,phi)
            if ~isempty(obj.ADL_up.Layers)
            obj.ADL_up.solve();
            end
            obj.f = f;
            obj.theta = theta;
            obj.phi = phi;
            Z_a = zeros(1,length(obj.f));
            
            mx = -20:20;
            my = -20:20;

            [MX, MY, F] = meshgrid(mx, my, obj.f);
                        
            K0 = 2 * pi * F ./ 3e8;
            KX = K0 * sin(obj.theta) * cos(obj.phi);
            KY = K0 * sin(obj.theta) * sin(obj.phi);

            KXM = KX - 2*pi.*MX ./ obj.dx;
            KYM = KY - 2*pi.*MY ./ obj.dy;

            KYM_DOWN = KYM;

            if obj.CavityWalls
                KYM_DOWN = - 2*pi.*MY ./ obj.dy;
            end

            KRO_UP = sqrt(KXM .^ 2 + KYM .^ 2);
            KRO_DOWN = sqrt(KXM.^2 + KYM_DOWN .^2);
            
            D_inf = obj.D_up(KXM, KYM, KRO_UP, F) + obj.D_down(KXM, KYM_DOWN, KRO_DOWN, F);
            Z_a = 1./obj.dx .* sum(-sinc(KXM(1,:,:) .* obj.delta_d ./ 2./ pi).^2./D_inf, 2);
            Z_a = (reshape(Z_a, [1,length(obj.f)]));
            obj.ADL_up.Frequencies = f;

        end
        
        function Z_a = activeInputImpedanceSynth(obj, f, theta ,phi, Z_up_TE, Z_up_TM)
            if ~isempty(obj.ADL_up.Layers)
%             obj.ADL_up.solve();
            end
            obj.f = f;
            obj.theta = theta;
            obj.phi = phi;
            
            mx = -10:10;
            my = -10:10;

            [MX, MY, F] = meshgrid(mx, my, obj.f);
                        
            K0 = 2 * pi * F ./ 3e8;
            KX = K0 * sin(obj.theta) * cos(obj.phi);
            KY = K0 * sin(obj.theta) * sin(obj.phi);

            KXM = KX - 2*pi.*MX ./ obj.dx;
            KYM = KY - 2*pi.*MY ./ obj.dy;

            KYM_DOWN = KYM;

            if obj.CavityWalls
                KYM_DOWN = - 2*pi.*MY ./ obj.dy;
            end

            KRO_UP = sqrt(KXM .^ 2 + KYM .^ 2);
            KRO_DOWN = sqrt(KXM.^2 + KYM_DOWN .^2);
            
            D_inf = obj.D_upSynth(KXM, KYM, KRO_UP, F, Z_up_TE, Z_up_TM) + obj.D_down(KXM, KYM_DOWN, KRO_DOWN, F);
            Z_a = 1./obj.dx .* sum(-sinc(KXM(1,:,:) .* obj.delta_d ./ 2./ pi).^2./D_inf, 2);
            Z_a = (reshape(Z_a, [1,length(obj.f)]));
            obj.ADL_up.Frequencies = f;

        end
        
        function [Z_up_TE, Z_up_TM] = getADLimpedance(obj, f, theta ,phi)
            if ~isempty(obj.ADL_up.Layers)
                obj.ADL_up.solve();
            end
            obj.f = f;
            obj.theta = theta;
            obj.phi = phi;
            
            mx = -10:10;
            my = -10:10;

            [MX, MY, F] = meshgrid(mx, my, obj.f);
                        
            K0 = 2 * pi * F ./ 3e8;
            KX = K0 * sin(obj.theta) * cos(obj.phi);
            KY = K0 * sin(obj.theta) * sin(obj.phi);

            KXM = KX - 2*pi.*MX ./ obj.dx;
            KYM = KY - 2*pi.*MY ./ obj.dy;
            KRO_UP = sqrt(KXM .^ 2 + KYM .^ 2);
            [Z_up_TE, Z_up_TM] = obj.ADL_up.getInputImpedance_ConnectedArray(F, obj.theta, obj.phi, KRO_UP);
            Z_up_TE(isnan(Z_up_TE)) = 0;
            Z_up_TM(isnan(Z_up_TM)) = 0;
        end
        
        
        
        function D = D_up(obj,KX,KY, KRO, F)
          
           [Z_up_TE, Z_up_TM] = obj.ADL_up.getInputImpedance_ConnectedArray(F, obj.theta, obj.phi, KRO);
%             [Z_up_TE, Z_up_TM] = obj.ADL_up.getInputImpedance(obj.f, obj.theta);
            Z_up_TE(isnan(Z_up_TE)) = 0;  %% is this the correct way to handle this?
            Z_up_TM(isnan(Z_up_TM)) = 0;  %% is this the correct way to handle this?


G_up = -1 .* ((1./Z_up_TE .* KX.^2 + 1./Z_up_TM .* KY.^2) ./ (KX.^2 + KY.^2));
            
            D = 1./obj.dy .* sum(G_up .* besselj(0, KY .* obj.w ./ 2), 1);
        end
        
        function D = D_upSynth(obj,KX,KY, ~, ~, Z_up_TE, Z_up_TM)
          
%             [Z_up_TE, Z_up_TM] = obj.ADL_up.getInputImpedance(obj.f, obj.theta);
            Z_up_TE(isnan(Z_up_TE)) = 0;  %% is this the correct way to handle this?
            Z_up_TM(isnan(Z_up_TM)) = 0;  %% is this the correct way to handle this?


            G_up = -1 .* ((1./Z_up_TE .* KX.^2 + 1./Z_up_TM .* KY.^2) ./ (KX.^2 + KY.^2));
            
            D = 1./obj.dy .* sum(G_up .* besselj(0, KY .* obj.w ./ 2), 1);
        end
        
        function D = D_down(obj, KX, KY,KRO, F)
             [Z_down_TE, Z_down_TM] = obj.Zin_down(KRO, F);
            G_down = -1 .* ((1./Z_down_TE .* KX.^2 + 1./Z_down_TM .* KY.^2) ./ (KX.^2 + KY.^2));
            
            D = 1./obj.dy .* sum(G_down .* besselj(0, KY .* obj.w ./ 2), 1);
        end
        
        function [Zin_TE, Zin_TM] = Zin_down(obj, KRO, F)
            % Transfer short through air section;
        
            % Air section
            k0 = 2 .* pi .* F ./ 3e8;  
            
            kz0 = -1j .* sqrt(- (k0.^2 - KRO.^2)); 

            Z_TE_air = 120 .* pi  .* k0 ./ kz0;
            Z_TM_air = 120 .* pi  .* kz0 ./ k0;
            
            % Dielectric section
            k_dielectric = sqrt(obj.er_back) .* 2 .* pi .* F ./ 3e8;  
            kz_dielectric = -1j .* sqrt(- (k_dielectric.^2 - KRO.^2)); 

            Z_TE_dielectric = 120 .* pi ./ sqrt(obj.er_back) .* k_dielectric ./ kz_dielectric;
            Z_TM_dielectric = 120 .* pi ./ sqrt(obj.er_back) .* kz_dielectric ./ k_dielectric;
            
            % Transfer short through air section
            if obj.h == inf
                Zin_air_TE = Z_TE_air;
                Zin_air_TM = Z_TM_air;
            else
            Zin_air_TE =  impedance_transfer_tx_line(Z_TE_air, 0, kz0, obj.h - obj.dielectric_thickness);
            Zin_air_TM =  impedance_transfer_tx_line(Z_TM_air, 0, kz0, obj.h - obj.dielectric_thickness);
            end
            % Transfer input impedance of air section through dielectric 
            Zin_TE =  impedance_transfer_tx_line(Z_TE_dielectric, Zin_air_TE, kz_dielectric, obj.dielectric_thickness);
            Zin_TM =  impedance_transfer_tx_line(Z_TM_dielectric, Zin_air_TM, kz_dielectric, obj.dielectric_thickness);
        end
    end
end

