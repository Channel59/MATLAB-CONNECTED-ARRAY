function [Z_in] = impedance_transfer_tx_line(Z0, Zl, k0, l)
%IMPEDANCE_TRANSFER_TX_LINE Summary of this function goes here
%   Detailed explanation goes here
Z_in = Z0 .* (Zl + 1j .* Z0 .* tan(k0 .* l)) ./ (Z0 + 1j .* Zl .* tan(k0 .* l));
end

