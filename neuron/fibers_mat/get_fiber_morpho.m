function fiber_morpho = get_fiber_morpho(fiber_model, diameter)
% Return the ultrastructural morphological features of a myelinated fiber model.
% Surrogate of part of the Python implementation build_fiber.py.
% 
%   fiber_morpho = get_fiber_morpho(fiber_model, diameter)
%   fiber_model     choices: "ascent_mrg"
%   diameter        diameter of the fiber (axon + myelin domains)
% 
%   ---------------------------------------------------------------
%   fiber_morpho (fields)       units   description
%   ---------------------------------------------------------------
%   diam                        [um]    diameter of fiber (axon+myelin)
%   diam_ranvier                [um]    diameter of Ranvier node
%   diam_mysa                   [um]    diameter of MYSA
%   diam_flut                   [um]    diameter of FLUT
%   diam_stin                   [um]    diameter of STIN (axon)
%   length_internode            [um]    length of the internode (including 1/2 of Ranvier nodes at each side)
%   length_ranvier              [um]    length of Ranvier node
%   length_mysa                 [um]    length of MYSA
%   length_flut                 [um]    length of FLUT
%   length_stin                 [um]    length of STIN
%   n_compartments_internode    [-]     number of compartments in each internode
%   mask_ranvier_internode      (bool)  mask encoding Ranvier node positions among the internode compartments
%   mask_mysa_internode         (bool)  mask encoding MYSA positions among the internode compartments
%   mask_flut_internode         (bool)  mask encoding FLUT positions among the internode compartments
%   mask_stin_internode         (bool)  mask encoding STIN positions among the internode compartments
%   n_lamellae                  [-]     number of lamellae in the myelin sheath
%   ---------------------------------------------------------------
% 
% See also get_fiber_pos_z

    switch fiber_model

        case 'ascent_mrg'
            
            length_mysa	= 3;
            length_ranvier = 1.0;
            
            g = [];
            diam_ranvier = 0.01093*diameter.^2 + 0.1008*diameter + 1.099;
            diam_mysa = diam_ranvier;
            diam_flut = 0.02361*diameter.^2 + 0.3673*diameter + 0.7122;
            diam_stin = diam_flut;
            
            length_internode = zeros(1,length(diameter));
            idcs_large_diameters = diameter >= 5.643;
            length_internode(idcs_large_diameters) = -8.215*diameter(idcs_large_diameters).^2 + 272.4*diameter(idcs_large_diameters) - 780.2;
            length_internode(~idcs_large_diameters) = 81.08*diameter(~idcs_large_diameters) + 37.84;
            
            length_flut = -0.1652*diameter.^2 + 6.354*diameter - 0.2862;
            n_lamellae = -0.4749*diameter.^2 + 16.85*diameter - 0.7648;
            
            length_stin = (length_internode - length_ranvier - (2*length_mysa) - (2*length_flut))/6;

            n_compartments_internode = 11;
            mask_ranvier_internode = [true false false false false false false false false false false];
            mask_mysa_internode = [false true false false false false false false false false true];
            mask_flut_internode = [false false true false false false false false false true false];
            mask_stin_internode = [false false false true true true true true true false false];

        otherwise
            error("fiber_model not valid")

    end

    fiber_morpho.diam = diameter;
    fiber_morpho.diam_ranvier = diam_ranvier;
    fiber_morpho.diam_mysa = diam_mysa;
    fiber_morpho.diam_flut = diam_flut;
    fiber_morpho.diam_stin = diam_stin;
    fiber_morpho.length_ranvier = length_ranvier*ones(1,length(diameter));
    fiber_morpho.length_mysa = length_mysa*ones(1,length(diameter));
    fiber_morpho.length_flut = length_flut;
    fiber_morpho.length_stin = length_stin;
    fiber_morpho.length_internode = length_internode;
    fiber_morpho.n_compartments_internode = n_compartments_internode;
    fiber_morpho.mask_ranvier_internode = mask_ranvier_internode;
    fiber_morpho.mask_mysa_internode = mask_mysa_internode;
    fiber_morpho.mask_flut_internode = mask_flut_internode;
    fiber_morpho.mask_stin_internode = mask_stin_internode;
    fiber_morpho.n_lamellae = n_lamellae;

end