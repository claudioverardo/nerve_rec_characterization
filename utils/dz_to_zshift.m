function [z_shift, z_shift_m] = dz_to_zshift(fiber_model, d_fiber, dz)
% Get spatial discretization of the internode of a myelinated fiber
    fiber_morpho = get_fiber_morpho(fiber_model, d_fiber);
    z_shift_m = [0:dz:fiber_morpho.length_internode/2 fiber_morpho.length_internode/2 flip(fiber_morpho.length_internode-(0:dz:fiber_morpho.length_internode/2))]'*1e-6; % m
    z_shift = z_shift_m/(fiber_morpho.length_internode*1e-6);
end