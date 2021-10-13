function D_inf = D_infty(k0, dy, wd, kxm, kym, type, hs)

 
    kro = sqrt(kxm .^ 2 + kym .^ 2);
    
switch type
    case 'strip'
        [v_TM, v_TE, i_TM, i_TE] = trxline_freespace(k0, kro, 0);
       
        [SGFej] = SpectralGFej(k0, 1, kxm, kym, v_TM, v_TE, i_TM, i_TE);
        SGF_xx = SGFej(:,:,1,1);
        J = besselj(0,kym .* wd ./ 2);
        accum = (SGF_xx .* J);
        D_inf = 1./dy .* sum(accum,1);
    case 'strip_reflector'
        fprintf("Case: strip_reflector\n")
        [v_TM, v_TE, i_TM, i_TE] = trxline_GroundSlab(k0, 1, hs ,kro, hs);
       
        [SGFej] = SpectralGFej(k0, 1, kxm, kym, v_TM, v_TE, i_TM, i_TE);
        SGF_xx = SGFej(:,:,1,1);
        J = besselj(0,kym .* wd ./ 2);
        accum = (SGF_xx .* J);
        D_inf = 1./dy .* sum(accum,1);
    case 'slot'
        fprintf("Case: slot\n")
        [~, ~, i_TM, i_TE] = trxline_freespace_magnetic(k0, kro, 0);
        SGF_xx = -(i_TE .* kxm .^ 2+ i_TM .* kym .^2) ./ (kxm .^ 2 + kym .^2);  %SGFhm_XX
        J = besselj(0,kym .* wd ./ 2);
        accum = (SGF_xx .* J);
        D_inf = 1./dy .* sum(accum,1);
    case 'slot_reflector'
        fprintf("Case: slot_reflector\n")
        [~, ~, i_TM, i_TE] = trxline_slot_groundplate(k0, kro, 0, hs);
        SGF_xx = -(i_TE .* kxm .^ 2+ i_TM .* kym .^2) ./ (kxm .^ 2 + kym .^2);  %SGFhm_XX
        J = besselj(0,kym .* wd ./ 2);
        accum = (SGF_xx .* J);
        D_inf = 1./dy .* sum(accum,1);
end


end