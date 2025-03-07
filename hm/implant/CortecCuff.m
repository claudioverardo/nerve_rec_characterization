classdef CortecCuff < Electrode
    % CortecCuff Subclass of Electrode representing a Cortec cuff electrode
    % (active sites with [partial-]ring shape).
    % 
    % CortecCuff Properties:
    %     id                        - (inherited from Electrode)
    % 
    %     --------------- shape parameters ----------------
    %     radius                    -
    %     radius_nerve              -
    %     length                    -
    %     thick                     -
    %     n_as                      -
    %     length_as                 -
    %     depth_as                  -
    %     width_as                  -
    %     pitch_as                  -
    %     
    %     ------------- insertion parameters --------------
    %     x                         -
    %     y                         -
    %     z                         -
    %     theta_z                   -
    %     
    %     ----------- encapsulation parameters ------------
    %     thick_encaps              -
    % 
    % CortecCuff Properties (Dependent):
    %     n_sites                   -
    %     alpha_as                  -
    %     z_as                      -
    %     phi_as                    -
    %     
    % CortecCuff Methods:
    %     CortecCuff                - constructor
    %     add_to_model              -
    %     assign_materials          -
    %     set_active_site_current   -
    %     get_selections            -

    properties
        % radius - radius of the electrode (double) [m]
        radius

        % radius_nerve - external radius of the nerve (double) [m]
        radius_nerve

        % length - length of the electrode (double) [m]
        length 

        % thick - thickness of the electrode (double) [m]
        thick

        % n_as - number of active sites (int)
        n_as

        % length_as - length of each active site (double) [m]
        length_as

        % depth_as - depth of each active site (double) [m]
        depth_as

        % width_as - width of each active site (double) [m]
        width_as

        % pitch_as - distance between adjacent active sites (double) [m]
        pitch_as
        
        % x - x-position of the cuff center (double) [m]
        x
        
        % y - y-position of the cuff center (double) [m]
        y
        
        % z - z-position of the cuff center (double) [m]
        z

        % theta_z - angle of insertion around z axis (double) [deg]
        theta_z

        % thick_encaps - thickness of the encapsulation tissue (double) [deg]
        % set to zero (ie, no encapsulation) by default
        thick_encaps
    end
    
    properties (Dependent)
        % n_sites - number of active sites (int)
        n_sites

        % alpha_as - angle of coverage of the nerve by each active site (double) [deg]
        alpha_as

        % z_as - z-position of the active sites (n_sites x 1 double) [m]
        % active sites are ordered first-to-last with increasing z-position
        z_as
        
        % phi_as - angular position (xy plane) of the active sites, at their center (n_sites x 1 double) [deg]
        phi_as
    end
    
    methods
        
        function NSites = get.n_sites(obj)
            % GET.N_SITES get the number of active sites
            %
            % RETURNS
            % -------
            % NSites (1x1 int)
            %     number of active sites

            NSites = obj.n_as;
        end
        
        function AlphaAs = get.alpha_as(obj)
            % GET.ALPHA_AS get the angle of coverage of the nerve (active sites)
            %
            % RETURNS
            % -------
            % AlphaAs (1x1 double)
            %     angle of coverage of the nerve (with active site)

            AlphaAs = obj.length_as / obj.radius / (2*pi) * 360;
        end
        
        function ZAs = get.z_as(obj)
            % GET.Z_AS get the z-position of the active sites
            %
            % RETURNS
            % -------
            % ZAs (n_asx1 double)
            %     z-positions of the active sites
            
            if obj.n_as == 1
                ZAs = 0;
            else
                ZAs = (-(obj.n_as-1)/2:(obj.n_as-1)/2)*obj.pitch_as;
            end
            ZAs = ZAs + obj.z;
        end

        function PhiAs = get.phi_as(obj)
            % GET.PHI_AS get the angular position of the active sites (at their center)
            % in the xy cross-section of the cuff
            %
            % RETURNS
            % -------
            % PhiAs (n_asx1 double)
            %     angular position of the active sites

            PhiAs = obj.alpha_as/2 + obj.theta_z;
        end

         function obj = CortecCuff(params)
            % CORTECCUFF instantiate the object with the parameters in
            % the structure params
            %
            % PARAMETERS
            % ----------
            % params (1x1 struct)
            %     the structure with the parameters
            %
            % RETURNS
            % -------
            % obj (1x1 CortecCuff)
            %    the modified caller object
            
            obj.id = params.id;
            obj.radius = params.radius;
            obj.radius_nerve = params.radius_nerve;
            obj.length = params.length;
            obj.thick = params.thick;
            obj.n_as = params.n_as;
            obj.length_as = params.length_as;
            obj.depth_as = params.depth_as;
            obj.width_as = params.width_as;
            if obj.n_as == 1
                params.pitch_as = [];
            end
            obj.pitch_as = params.pitch_as;

            if ~isfield(params, 'x') || isempty(params.x) % legacy
                params.x = 0;
            end
            obj.x = params.x;
            if ~isfield(params, 'y') || isempty(params.y) % legacy
                params.y = 0;
            end
            obj.y = params.y;
            obj.z = params.z;
            obj.theta_z = params.theta_z;

            if ~isfield(params, 'thick_encaps') || isempty(params.thick_encaps)
                params.thick_encaps = 0;
            end
            obj.thick_encaps = params.thick_encaps;
            if obj.thick_encaps > 0
                warning('Rototranslation of electrodes is not currently supported with encapsulation');
            end
         end
        
        function model = add_to_model(obj, model)
            % ADD_TO_MODEL add the electrode to the model
            %
            % PARAMETERS
            % ----------
            % model (1x1 COMSOL ModelClient)
            %     the model to modify
            %
            % RETURNS
            % -------
            % model (1x1 COMSOL ModelClient)
            %     the modified COMSOL model

            geom1 = model.geom('geom1');
            
            % ======================== CREATE ELECTRODE BODY ========================
            elec_out = geom1.feature.create(strcat('elec_out_', int2str(obj.id)), 'Cylinder');
            elec_out.set('h', obj.length);
            elec_out.set('r', obj.radius + obj.thick);
            elec_out.set('pos', [0, 0, -obj.length/2]);
            
            % geom1.run()
            
            elec_in = geom1.feature.create(strcat('elec_in_', int2str(obj.id)), 'Cylinder');
            elec_in.set('h', obj.length);
            elec_in.set('r', obj.radius);
            elec_in.set('pos', [0, 0, -obj.length/2]);
            
            % geom1.run();

            elec_full = geom1.create(strcat('elec_full_', int2str(obj.id)), 'Difference');
            elec_full.selection('input').set(strcat('elec_out_', int2str(obj.id)));
            elec_full.selection('input2').set(strcat('elec_in_', int2str(obj.id)));
            
            % ========================= CREATE ACTIVE SITES =========================
            for j=1:obj.n_as

                as_coverage_wp = geom1.feature.create(strcat('as_coverage_wp_', int2str(obj.id), '_', int2str(j)), 'WorkPlane');
                as_coverage_wp.set('planetype', 'quick');
                as_coverage_wp.set('quickplane', 'xy');
                as_coverage_wp.set('quickz', -obj.length/2);
                
                % create mask used to draw the active site
                as_coverage_xy = as_coverage_wp.geom.feature.create(strcat("as_coverage_xy_", int2str(obj.id), '_', int2str(j)), 'Circle');
                as_coverage_xy.set('r', obj.radius + obj.thick);
                as_coverage_xy.set('angle', obj.alpha_as);
                % as_coverage_xy.set('rot', (2*pi-obj.alpha_as) / 2);
                as_coverage_xy.set('rot', 0);
                
                geom1.run(strcat('as_coverage_wp_', int2str(obj.id), '_', int2str(j)));
                
                as_coverage = geom1.feature.create(strcat('as_coverage_', int2str(obj.id), '_', int2str(j)), 'Extrude');
                as_coverage.selection('input').init();
                as_coverage.selection('input').set(strcat('as_coverage_wp_', int2str(obj.id), '_', int2str(j)));
                as_coverage.set('distance', obj.length);

                % geom1.run();

                as_out = geom1.feature.create(strcat('as_out_', int2str(obj.id), '_', int2str(j)), 'Cylinder');
                as_out.set('h', obj.width_as);
                as_out.set('r', obj.radius + obj.depth_as);
                as_out.set('pos', [0, 0, -obj.width_as/2 + obj.z_as(j) - obj.z]);
    
                % geom1.run();
    
                as_in = geom1.feature.create(strcat('as_in_', int2str(obj.id), '_', int2str(j)), 'Cylinder');
                as_in.set('h', obj.width_as);
                as_in.set('r', obj.radius);
                as_in.set('pos', [0, 0, -obj.width_as/2 + obj.z_as(j) - obj.z]);
    
                % geom1.run();
    
                as_round = geom1.create(strcat('as_round_', int2str(obj.id), '_', int2str(j)), 'Difference');
                as_round.selection('input').set(strcat('as_out_', int2str(obj.id), '_', int2str(j)));
                as_round.selection('input2').set(strcat('as_in_', int2str(obj.id), '_', int2str(j)));
    
                as_origin(j) = geom1.create(strcat('as_origin_', int2str(obj.id), '_', int2str(j)), 'Intersection');
                as_origin(j).selection('input').set({strcat('as_round_', int2str(obj.id), '_', int2str(j)), strcat('as_coverage_', int2str(obj.id), '_', int2str(j))});
    
                % geom1.run()

            end
            
            % =============== APPLY ROTO-TRANSLATION TO ELECTRODE BODY ==============
            elec_moved = geom1.create(strcat('elec_moved_', int2str(obj.id)), 'Move');
            elec_moved.selection('input').set(strcat('elec_full_', int2str(obj.id)));
            elec_moved.set('displ', [obj.x, obj.y, obj.z]);

            elec_rotz = geom1.feature.create(strcat('elec_rotz_', int2str(obj.id)), 'Rotate');
            elec_rotz.selection('input').set(strcat('elec_moved_', int2str(obj.id)));
            elec_rotz.set('pos', [0, 0, 0]);
            elec_rotz.set('axis', [0, 0, 1]);
            elec_rotz.set('rot', obj.theta_z);

            % =============== APPLY ROTO-TRANSLATION TO ACTIVE SITES ================
            as_names = {};
            for j=1:obj.n_as
            
                as_moved = geom1.create(strcat('as_moved_', int2str(obj.id), '_', int2str(j)), 'Move');
                as_moved.selection('input').set(strcat('as_origin_', int2str(obj.id), '_', int2str(j)));
                as_moved.set('displ', [obj.x, obj.y, obj.z]);

                as(j) = geom1.feature.create(strcat('as_', int2str(obj.id), '_', int2str(j)), 'Rotate');
                as(j).selection('input').set(strcat('as_moved_', int2str(obj.id), '_', int2str(j)));
                as(j).set('pos', [0, 0, 0]);
                as(j).set('axis', [0, 0, 1]);
                as(j).set('rot', obj.theta_z);
                as(j).set('createselection', 'on');
    
                as_names = [as_names strcat('as_', int2str(obj.id), '_', int2str(j))];

            end

            % ============== SUBTRACT ACTIVE SITES FROM ELECTRODE BODY ==============
            if ~isempty(as_names)
                elec = geom1.create(strcat('elec_', int2str(obj.id)), 'Difference');
                elec.selection('input').set(strcat('elec_rotz_', int2str(obj.id)));
                elec.selection('input2').set(as_names);
                elec.set("keepsubtract", true);
                elec.set('createselection', 'on');
            else
                elec = geom1.create(strcat('elec_', int2str(obj.id)), 'Copy');
                elec.selection('input').set(strcat('elec_rotz_', int2str(obj.id)));
                elec.set('createselection', 'on');
            end

            % ===================== CREATE ENCAPSULATION TISSUE =====================
            if obj.thick_encaps>0
                encaps_out = geom1.feature.create(strcat('encaps_out_', int2str(obj.id)), 'Cylinder');
                encaps_out.set('h', obj.length + 2*obj.thick_encaps);
                encaps_out.set('r', obj.radius + obj.thick + obj.thick_encaps);
                encaps_out.set('pos', [0, 0, -(obj.length+2*obj.thick_encaps)/2]);
                
                % geom1.run()
                
                encaps_in = geom1.feature.create(strcat('encaps_in_', int2str(obj.id)), 'Cylinder');
                encaps_in.set('h', obj.length + 2*obj.thick_encaps);
                encaps_in.set('r', obj.radius - obj.thick_encaps);
                encaps_in.set('pos', [0, 0, -(obj.length+2*obj.thick_encaps)/2]);
                
                % geom1.run();
    
                encaps_full = geom1.create(strcat('encaps_full_', int2str(obj.id)), 'Difference');
                encaps_full.selection('input').set(strcat('encaps_out_', int2str(obj.id)));
                encaps_full.selection('input2').set(strcat('encaps_in_', int2str(obj.id)));
                
                % geom1.run();
    
                encaps = geom1.create(strcat('encaps_', int2str(obj.id)), 'Difference');
                encaps.selection('input').set(encaps_full.tag.char);
                encaps.selection('input2').set([{elec.tag.char} as_names]);
                encaps.set("keepsubtract", true);
                encaps.set('createselection', 'on');
            end
            
            geom1.run()
            
            % =============== CREATE CURRENT TERMINALS IN THE PHYSICS ===============
            for j=1:obj.n_as
                assel = mphgetselection(model.selection(strcat('geom1_as_', int2str(obj.id), '_', int2str(j), '_dom')));
                asdom = assel.entities;
                model.physics('ec').create(strcat('term_', int2str(obj.id), '_', int2str(j)), 'DomainTerminal', 3);
                model.physics('ec').feature(strcat('term_', int2str(obj.id), '_', int2str(j))).label(strcat('term_', int2str(obj.id), '_', int2str(j)));
                model.physics('ec').feature(strcat('term_', int2str(obj.id), '_', int2str(j))).selection.set(asdom);
            end
        end

        function model = assign_materials(obj, model)
            % ASSIGN_MATERIALS assign materials to electrode shaft/active sites
            %
            % PARAMETERS
            % ----------
            % model (1x1 COMSOL ModelClient)
            %     the model to modify
            %
            % RETURNS
            % -------
            % model (1x1 COMSOL ModelClient)
            %     the modified COMSOL model

            elecdom = mphgetselection(model.material('elec').selection()).entities;
            asdom = mphgetselection(model.material('as').selection()).entities;

            elecsel = mphgetselection(model.selection(strcat('geom1_elec_', int2str(obj.id), '_dom')));
            elecdom = [elecdom, elecsel.entities];

            for k=1:obj.n_as
                assel = mphgetselection(model.selection(strcat('geom1_as_', int2str(obj.id), '_', int2str(k), '_dom')));
                asdom = [asdom, assel.entities];
            end

            model.material('elec').selection.set(elecdom);
            model.material('as').selection.set(asdom);

            if obj.thick_encaps>0
                encapsdom = mphgetselection(model.material('encapsulation').selection()).entities;
                encapssel = mphgetselection(model.selection(strcat('geom1_encaps_', int2str(obj.id), '_dom')));
                encapsdom = [encapsdom, encapssel.entities];
                model.material('encapsulation').selection.set(encapsdom);
            end
        end
        
        function model = set_active_site_current(obj, model, site_id, curr)
            % SET_ACTIVE_SITE_CURRENT set the current of an active site
            %
            % PARAMETERS
            % ----------
            % model (1x1 COMSOL ModelClient)
            %     the model to modify
            % site_id (1x1 integer)
            %     active site numerical identifier
            % curr (1x1 double)
            %     current amplitude
            %
            % RETURNS
            % -------
            % model (1x1 COMSOL ModelClient)
            %     the modified COMSOL model
            
            model.physics('ec').feature(strcat('term_', int2str(obj.id), '_', int2str(site_id))).set('I0', curr);
        end
        
        function [elec_sel, as_sel, encaps_sel] = get_selections(obj)
            elec_sel = {strcat('elec_', num2str(obj.id))};
            as_sel = {};
            for k = 1:obj.n_as
                as_sel{end+1} = strcat('as_', int2str(obj.id), '_', int2str(k));
            end
            if obj.thick_encaps>0
                encaps_sel = {strcat('encaps_', num2str(obj.id))};
            else
                encaps_sel = {};
            end
        end
        
    end
end