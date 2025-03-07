function [fiber_pos_z, fiber_mask_ranvier, fiber_morpho] = get_fiber_pos_z(fiber_params)
% Build the z-positions of nodes (compartments) in a straight fiber. 
% Surrogate of part of the Python implementation build_fiber.py.
% 
%   [fiber_pos_z, fiber_mask_ranvier, fiber_morpho] = get_fiber_pos_z(fiber_params)
%   fiber_pos_z             [mm]
%   fiber_mask_ranvier      bool
%   fiber_morpho            (struct) see get_fiber_morpho()
%   fiber_params            (struct) see below
% 
%   ---------------------------------------------------------------
%   fiber_params (fields)   myel    unmyel  description
%   ---------------------------------------------------------------
%   fiber_model*            x       x       (string)
%   diameter                x       x       (double) [um]
%   n_internodes            x               (int)
%   n_nodes                         x       (int)
%   node_length                     x       (double)
%   center_node_location    x       x       (3x1 double) [mm]
%   ---------------------------------------------------------------
% 
%   (*) see get_fiber_model_id() for the available fiber models

    if ~isfield(fiber_params, 'center_node_location')
        fiber_params.center_node_location = [0,0,0];
    end

    if check_fiber_myelinated(fiber_params.fiber_model)
        fiber_morpho = get_fiber_morpho(fiber_params.fiber_model, fiber_params.diameter);
        fiber_morpho.n_internodes = fiber_params.n_internodes;
              
        node_spacing_internode = [ % um
            (fiber_morpho.length_ranvier + fiber_morpho.length_mysa) / 2
            (fiber_morpho.length_mysa + fiber_morpho.length_flut) / 2
            (fiber_morpho.length_flut + fiber_morpho.length_stin) / 2
            fiber_morpho.length_stin
            fiber_morpho.length_stin
            fiber_morpho.length_stin
            fiber_morpho.length_stin
            fiber_morpho.length_stin
            (fiber_morpho.length_stin + fiber_morpho.length_flut) / 2
            (fiber_morpho.length_flut + fiber_morpho.length_mysa) / 2
            (fiber_morpho.length_mysa + fiber_morpho.length_ranvier) / 2
        ];
        node_spacing_fiber = repmat(node_spacing_internode, fiber_params.n_internodes, 1);
        fiber_pos_z = cumsum(node_spacing_fiber);
        fiber_pos_z = [0; fiber_pos_z];
        fiber_pos_z = fiber_pos_z - fiber_morpho.length_internode * (fiber_params.n_internodes/2);
        fiber_pos_z = fiber_pos_z / 1000 + fiber_params.center_node_location(3); % mm

        fiber_mask_ranvier = [repmat(fiber_morpho.mask_ranvier_internode', fiber_params.n_internodes, 1); true];

    elseif check_fiber_unmyelinated(fiber_params.fiber_model)
        fiber_morpho = [];
        length = fiber_params.n_nodes * fiber_params.node_length/1000; % mm

        nodes_spacing = fiber_params.node_length * ones(fiber_params.n_nodes-1,1);
        fiber_pos_z = cumsum(nodes_spacing)/1000;
        fiber_pos_z = [0; fiber_pos_z];
        fiber_pos_z = fiber_pos_z - (length - fiber_params.node_length/1000)/2;
        fiber_pos_z = fiber_pos_z + fiber_params.center_node_location(3); % mm

        fiber_mask_ranvier = false(size(fiber_pos_z));

    else
        error("not valid fiber_model");
    end

end