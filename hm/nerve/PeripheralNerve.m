classdef PeripheralNerve
    % PeripheralNerve Nerve w/ fascicles (epineurium + perineurium + endoneurium).
    % 
    % PeripheralNerve Properties:
    %     length                            - length of the nerve
    %     radius                            - radius of the nerve epineurium
    %     fascicle_topography               - fascicles in the nerve
    %     fiber_topography                  - fibers in the nerve
    % 
    % PeripheralNerve Methods:
    %     PeripheralNerve                   - constructor
    %     add_to_model                      -
    %     assign_materials                  -
    %     get_node_locs                     -

    properties
        % length - length of the nerve (double) [m]
        length
        
        % radius - radius of the nerve epineurium (double) [m]
        radius
                                    
        % fascicle_topography - fascicular topography of the nerve (1 x 1 FascicleTopography)
        fascicle_topography
        
        % fiber_topography - functional topography of the nerve (1 x 1 FiberTopography)
        fiber_topography
    end

    methods

        function obj = PeripheralNerve(params)
            % Instantiate object properties from struct of parameters
            % 
            % obj = PeripheralNerve(params)

            obj.length = params.nerve_extrusion_length;
            obj.radius = params.radius_nerve;
            
            % fascicular topography (init object)
            obj.fascicle_topography = FascicleTopography(params);

            % fiber topography (init object)
            params.fibers.dir_model = params.dir_model;
            obj.fiber_topography = FiberTopography(params.fibers);
        end
        
        function model = add_to_model(obj, model)
            % model = add_to_model(model)

            geom1 = model.geom('geom1');

            % create epineurium
            disp('Adding nerve boundary')

            x = 0;
            y = 0;
            r = obj.radius;
            
            epifull = geom1.feature.create('epifull', 'Cylinder');
            epifull.set('h', obj.length);
            epifull.set('x', x);
            epifull.set('y', y);
            epifull.set('r', r);
            epifull.set('z', -obj.length/2);  % centers in z = 0

            geom1.run;
    
            epicopy = geom1.feature.create('epicopy', 'Copy');
            epicopy.selection('input').set('epifull');

            salgeom = geom1.feature.create('salgeom','Difference');
            salgeom.selection('input').set('salfull');
            salgeom.selection('input2').set('epicopy');
            salgeom.set('createselection','on');

            disp('Adding fascicles')
            model = obj.fascicle_topography.add_to_model(model);
        end
        
        function model = assign_materials(obj, model)
            % model = assign_materials(model)

            episel = mphgetselection(model.selection('geom1_epigeom_dom'));
            epidom = episel.entities;
            model.material('epi').selection.set(epidom);
            
            endosel = mphgetselection(model.selection('geom1_endogeom_dom'));
            endodom = endosel.entities;
            model.material('endo').selection.set(endodom);

        end

        function [node_locs, node_locs_metadata] = get_node_locs(obj, id_fibers, varargin)
            % Generate 3D locations of fiber nodes (compartments) from fibers in the 2D functional topography.
            % NOTE: wrapper to FiberTopography.get_node_locs() to avoid specifying the nerve length.
            % 
            % [node_locs, node_locs_metadata] = get_node_locs(id_fibers, varargin)
            % 
            % PARAMETERS
            % ----------
            % id_fibers (n_fibers x 1 int)
            %     ids of the fibers 
            % 
            % NAME, VALUE ARGS
            % ----------------
            % node_length_unmyel (double)
            %     spatial discretization step for unmyelinated fibers (50/6 um by default)
            % log_step (int)
            %     control amount of log information in the command window
            %
            % RETURNS
            % -------
            % node_locs (n_fibers x 3 double)
            %     content of locs.txt files (see HMModel.export_node_locs)
            % node_locs_metadata (struct)
            %     content of locs_meta.mat files (see HMModel.export_node_locs)

            if ~exist("id_fibers",'var') || isempty(id_fibers)
                id_fibers = obj.fiber_topography.topo_func(:,1);
            end

            [node_locs, node_locs_metadata] = obj.fiber_topography.get_node_locs(id_fibers, obj.length, varargin{:});
        end
        
    end
end