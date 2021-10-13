classdef ConnectedArray < handle
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
        function obj = ConnectedArray(dx, dy, h, w, delta_d,er, dielectric_thickness)
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
            obj.theta = theta;
            obj.phi = phi;
            Z_a = zeros(1,length(f));
            for i = 1:length(f)
                obj.f = f(i);
                obj.ADL_up.Frequencies = f(i);
                %METHOD1 Summary of this method goes here
                %   Detailed explanation goes here
                mx = -20:20;
                my = -20:20;

                k = 2 * pi * obj.f ./ 3e8;

                k_cart = k_to_cartesian(k, obj.theta, obj.phi);
                kxm = k_cart(:,:,1) - 2 .* pi .* mx ./ obj.dx;
                kym = k_cart(:,:,2) - 2 .* pi .* my ./ obj.dy;

                if ~obj.CavityWalls
                    kym_down = kym;
                else
                    kym_down = - 2 * pi * my ./ obj.dy;
                end
                [KXM, KYM] = meshgrid(kxm, kym);
                [~, KYM_DOWN] = meshgrid(kxm, kym_down);

                D_inf = obj.D_up(KXM, KYM) + obj.D_down(KXM, KYM_DOWN);
                Z_a(i) = 1./obj.dx .* sum(-sinc(kxm .* obj.delta_d ./ 2./ pi).^2./D_inf);
            end
        end
        
        function D = D_up(obj, KX, KY)
            KRO = sqrt(KX .^ 2 + KY .^ 2);
            

            [Z_up_TE, Z_up_TM] = obj.ADL_up.getInputImpedance_ConnectedArray(obj.f, obj.theta, obj.phi, KRO);
%             [Z_up_TE, Z_up_TM] = obj.ADL_up.getInputImpedance(obj.f, obj.theta);
            Z_up_TE(isnan(Z_up_TE)) = 0;  %% is this the correct way to handle this?
            Z_up_TM(isnan(Z_up_TM)) = 0;  %% is this the correct way to handle this?


G_up = -1 .* ((1./Z_up_TE .* KX.^2 + 1./Z_up_TM .* KY.^2) ./ (KX.^2 + KY.^2));
            
            D = 1./obj.dy .* sum(G_up .* besselj(0, KY .* obj.w ./ 2), 1);
        end
        
        function D = D_down(obj, KX, KY)
            KRO = sqrt(KX .^ 2 + KY .^ 2);
             [Z_down_TE, Z_down_TM] = obj.Zin_down(KRO);
            G_down = -1 .* ((1./Z_down_TE .* KX.^2 + 1./Z_down_TM .* KY.^2) ./ (KX.^2 + KY.^2));
            
            D = 1./obj.dy .* sum(G_down .* besselj(0, KY .* obj.w ./ 2), 1);
        end
        
        function [Zin_TE, Zin_TM] = Zin_down(obj, KRO)
            % Transfer short through air section;
        
            % Air section
            k0 = 2 .* pi .* obj.f ./ 3e8;  
            
            kz0 = -1j .* sqrt(- (k0.^2 - KRO.^2)); 

            Z_TE_air = 120 .* pi  .* k0 ./ kz0;
            Z_TM_air = 120 .* pi  .* kz0 ./ k0;
            
            % Dielectric section
            k_dielectric = sqrt(obj.er_back) .* 2 .* pi .* obj.f ./ 3e8;  
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

