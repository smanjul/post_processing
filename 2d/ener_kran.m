function [E] = ener_kran(w_amp,mod_no_max,n1)

% Sends back contribution to KE by w_amp
% w_amp: Velocity component transformed to spectal space for one particular
% z. Has dimension 1 x n2. 
% mod_no_max: maximum azimuthal wavenumbers for which E is calculated
E = zeros(mod_no_max+1,1);

for l = 1:mod_no_max
    k_pos = l + 1;
    k_neg = n1 - l;
    amp_pk = w_amp(k_pos,:); % amplitude corresponding to k. 
    amp_nk = w_amp(k_neg,:); % amplitude corresponding to -k.
    E(l+1) = 0.5*(amp_nk*amp_nk' + amp_pk*amp_pk');
end

% For mean
k_pos = 1;
amp_pk = w_amp(k_pos,:); 
E(1) = 0.5*amp_pk*amp_pk';