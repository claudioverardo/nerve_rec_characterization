classdef uNG < Electrode
    % uNG Subclass of Electrode representing a uNG needle electrode
    % (uNG: microneurography).
    % 
    % uNG Properties:
    %     id                        - (inherited from Electrode)
    % 
    %     ----------- electrode shaft parameters ----------
    %     h_shaft                   -
    %     d_shaft                   -
    % 
    %     ------------ electrode tip parameters -----------
    %     slope_tip                 -
    %     d_rounding_tip            -
    %     h_uncoated_tip            -
    % 
    %     --------- insertion location parameters ---------
    %     x                         -
    %     y                         -
    %     z                         -
    % 
    %     -------- insertion orientation parameters -------
    %     theta_x                   -
    %     theta_y                   -
    %     theta_z                   -
    % 
    % uNG Properties (Dependent):
    %     n_sites                   -
    %     delta_rounding_tip        -
    %     h_tip                     -
    %     d_uncoated_tip            -
    %     
    % uNG Methods:
    %     uNG                       - constructor
    %     add_to_model              -
    %     assign_materials          -
    %     set_active_site_current   -
    %     get_selections            -

    properties
        % h_shaft - height of the electrode shaft (double) [m]
        h_shaft

        % d_shaft - diameter of the electrode shaft (double) [m]
        d_shaft

        % slope_tip - slope of the tapering of the electrode tip (double) [-]
        % slope_tip = (h_tip + delta_rounding_tip) / (d_shaft/2)
        slope_tip

        % d_rounding_tip - diameter of the rounding at the electrode tip (double) [m]
        d_rounding_tip

        % h_uncoated_tip - height of the uncoated part of the electrode tip (double) [m]
        h_uncoated_tip

        % x - origin x-location of the electrode reference frame (position of the tip vertex) (double)  [m]
        x

        % y - origin y-location of the electrode reference frame (position of the tip vertex) (double) [m]
        y

        % z - origin z-location of the electrode reference frame (position of the tip vertex) (double) [m]
        z
        
        % theta_x - orientation around x-axis of the electrode reference frame (double) [deg]
        theta_x
        
        % theta_y - orientation around y-axis of the electrode reference frame (double) [deg]
        theta_y
        
        % theta_z - orientation around z-axis of the electrode reference frame (double) [deg]
        theta_z
    end
    
    properties (Dependent)
        % n_sites - number of active sites (hardcoded to 1) (int)
        n_sites

        % delta_rounding_tip - height differenct wrt the electrode w/o rounding at the tip (double) [m]
        delta_rounding_tip

        % h_tip - height of the tip (uncoated+coated parts) (double) [m]
        h_tip

        % d_uncoated_tip - diameter of the uncoated part of the electrode tip (double) [m]
        d_uncoated_tip
    end
    
    methods
        
        function NSites = get.n_sites(obj)
            NSites = 1;
        end

        function DeltaRoundingTip = get.delta_rounding_tip(obj)
            ThetaRoundingTip = atan(1/obj.slope_tip);
            DeltaRoundingTip = obj.d_rounding_tip/2 * (sin(ThetaRoundingTip) + obj.slope_tip*cos(ThetaRoundingTip) - 1);
        end

        function HTip = get.h_tip(obj)
            HTip = obj.d_shaft/2 * obj.slope_tip - obj.delta_rounding_tip;
        end

        function DUncoatedTip = get.d_uncoated_tip(obj)
            DUncoatedTip = 2*(obj.h_uncoated_tip+obj.delta_rounding_tip)/obj.slope_tip;
        end
        
        function obj = uNG(params)
            obj.id = params.id;
            obj.h_shaft = params.h_shaft;
            obj.d_shaft = params.d_shaft;
            obj.slope_tip = params.slope_tip;
            obj.d_rounding_tip = params.d_rounding_tip;
            obj.h_uncoated_tip = params.h_uncoated_tip;
            obj.x = params.x;
            obj.y = params.y;
            obj.z = params.z;
            obj.theta_x = params.theta_x;
            obj.theta_y = params.theta_y;
            obj.theta_z = params.theta_z;
        end
        
        function model = add_to_model(obj, model) 
            geom1 = model.geom('geom1');

            %% Generate electrode body in electrode reference frame
            disp('Generating electrode geometry in electrode reference frame')

            % create domain corresponding to the shaft of the electrode
            elec_shaft = geom1.create(strcat('elec_shaft_', int2str(obj.id)), 'Cylinder');
            % elec_shaft.label('Cylinder uNG shaft');
            elec_shaft.set('pos', [0 0-obj.h_shaft-obj.h_tip 0]);
            elec_shaft.set('axis', [0 1 0]);
            elec_shaft.set('r', obj.d_shaft/2);
            elec_shaft.set('h', obj.h_shaft);
            elec_shaft.set('createselection', 'on');

            elec_tip_workplane = geom1.create(strcat('elec_tip_workplane_', int2str(obj.id)), 'WorkPlane');
            elec_tip_workplane.set('quickz', 0);
            elec_tip_workplane.set('unite', true);

            elec_tip_workplane_pol = elec_tip_workplane.geom.create(strcat(elec_tip_workplane.tag.char, '_pol'), 'Polygon');
            % elec_tip_workplane_pol.label('Polygon uNG_tip cross-section');
            elec_tip_workplane_pol.set('source', 'table');
            elec_tip_workplane_pol.set('table', ...
                [0 - obj.d_shaft/2      0 - obj.h_tip; ...
                0 + obj.d_shaft/2       0 - obj.h_tip; ...
                0                       0 + obj.delta_rounding_tip; ...
                0 - obj.d_shaft/2       0 - obj.h_tip]);

            elec_tip_workplane_fil = elec_tip_workplane.geom.create(strcat(elec_tip_workplane.tag.char, '_fil'), 'Fillet');
            % elec_tip_workplane_fil.label('Fillet uNG_tip');
            elec_tip_workplane_fil.set('radius', obj.d_rounding_tip/2);
            elec_tip_workplane_fil.selection('point').set(elec_tip_workplane_pol.tag.char, 2);
            
            elec_tip_workplane_r = elec_tip_workplane.geom.create(strcat(elec_tip_workplane.tag.char, '_r'), 'Rectangle');
            % elec_tip_workplane_r.label('Rectangle mask half uNG_tip section');
            elec_tip_workplane_r.set('pos', [0 0-obj.h_tip]);
            elec_tip_workplane_r.set('size', [obj.d_shaft/2 obj.h_tip]);

            elec_tip_workplane_int = elec_tip_workplane.geom.create(strcat(elec_tip_workplane.tag.char, '_int'), 'Intersection');
            elec_tip_workplane_int.selection('input').set({elec_tip_workplane_fil.tag.char elec_tip_workplane_r.tag.char});

            elec_tip_rev = geom1.create(strcat('elec_tip_rev_', int2str(obj.id)), 'Revolve');
            % elec_tip_rev.label('Revolve uNG tip cross-section');
            elec_tip_rev.set('angtype', 'full');
            elec_tip_rev.set('pos', [0 0]);
            elec_tip_rev.selection('input').set({elec_tip_workplane.tag.char});
            elec_tip_rev.set('createselection', 'on');

            %% Generate active site (uncoated tip) in electrode reference frame

            elec_tip_uncoated_block = geom1.create(strcat('elec_tip_uncoated_block_', int2str(obj.id)), 'Block');
            % elec_tip_uncoated_block.label('Block mask uncoated tip');
            elec_tip_uncoated_block.set('pos', [0 0-obj.h_uncoated_tip/2 0]);
            elec_tip_uncoated_block.set('axis', [0 1 0]);
            elec_tip_uncoated_block.set('base', 'center');
            elec_tip_uncoated_block.set('size', [obj.d_uncoated_tip obj.d_uncoated_tip obj.h_uncoated_tip]);

            % create domain corresponding to the uncoated tip of the electrode
            % tag 'elec_as_#' to be consistent w/ other implants
            elec_tip_uncoated = geom1.create(strcat('as_pre_Rt_', int2str(obj.id)), 'Intersection');
            % elec_tip_uncoated.label('Intersection mask uncoated tip - tip');
            elec_tip_uncoated.set('keep', true);
            elec_tip_uncoated.selection('input').set({elec_tip_uncoated_block.tag.char elec_tip_rev.tag.char});
            elec_tip_uncoated.set('createselection', 'on');

            elec_tip_uncoated_del = geom1.create(strcat('elec_tip_uncoated_del_', int2str(obj.id)), 'Delete');
            % elec_tip_uncoated_del.label('Delete block mask uncoated tip');
            elec_tip_uncoated_del.selection('input').init(3); 
            elec_tip_uncoated_del.selection("input").set(elec_tip_uncoated_block.tag.char, 1);

            % create domain corresponding to the coated tip of the electrode
            elec_tip_coated = geom1.create(strcat('elec_tip_coated_', int2str(obj.id)), "Difference");
            % elec_tip_coated.label('Difference tip - uncoated tip');
            elec_tip_coated.selection("input").set(elec_tip_rev.tag.char);
            elec_tip_coated.selection("input2").set(elec_tip_uncoated.tag.char);
            elec_tip_coated.set("keepsubtract", true);
            elec_tip_coated.set("keepadd", false);
            elec_tip_coated.set('createselection', 'on');

            % create domain corresponding to the coated part of the electrode (shaft+tip)
            % tag 'elec_#' to be consistent w/ other implants
            elec_pre_Rt = geom1.create(strcat('elec_pre_Rt_', int2str(obj.id)), "Union");
            % elec_pre_Rt.label('Union coated tip - shaft');
            elec_pre_Rt.selection("input").set({elec_tip_coated.tag.char, elec_shaft.tag.char});
            elec_pre_Rt.set("keep", false);
            elec_pre_Rt.set('createselection', 'on');

            %% Moving electrode geometry to laboratory reference frame
            disp('Moving electrode geometry to laboratory reference frame')

            % Apply roto-translation to electrode body
            rotz = geom1.feature.create(strcat('elec_rotz_', int2str(obj.id)),'Rotate');
            rotz.selection('input').set(elec_pre_Rt.tag.char);
            rotz.set('pos', [0, 0, 0]);
            rotz.set('axis', [0, 0, 1]);
            rotz.set('rot', obj.theta_z);
            
            roty = geom1.feature.create(strcat('elec_roty_', int2str(obj.id)),'Rotate');
            roty.selection('input').set(strcat('elec_rotz_', int2str(obj.id)));
            roty.set('pos', [0, 0, 0]);
            roty.set('axis', [0, 1, 0]);
            roty.set('rot', obj.theta_y);
            
            rotx = geom1.feature.create(strcat('elec_rotx_', int2str(obj.id)),'Rotate');
            rotx.selection('input').set(strcat('elec_roty_', int2str(obj.id)));
            rotx.set('pos', [0, 0, 0]);
            rotx.set('axis', [1, 0, 0]);
            rotx.set('rot',obj.theta_x);
            
            v = [obj.x, obj.y, obj.z];
            elec = geom1.feature.create(strcat('elec_', int2str(obj.id)),'Move');
            elec.selection('input').set(strcat('elec_rotx_', int2str(obj.id)));
            elec.set('displ', v);
            elec.set('createselection', 'on');
            
            % Apply roto-translation to electrode active sites
            rotz = geom1.feature.create(strcat('as_rotz_', int2str(obj.id)),'Rotate');
            rotz.selection('input').set(elec_tip_uncoated.tag.char);
            rotz.set('pos', [0, 0, 0]);
            rotz.set('axis', [0, 0, 1]);
            rotz.set('rot', obj.theta_z);

            roty = geom1.feature.create(strcat('as_roty_', int2str(obj.id)),'Rotate');
            roty.selection('input').set(strcat('as_rotz_', int2str(obj.id)));
            roty.set('pos', [0, 0, 0]);
            roty.set('axis', [0, 1, 0]);
            roty.set('rot', obj.theta_y);

            rotx = geom1.feature.create(strcat('as_rotx_', int2str(obj.id)),'Rotate');
            rotx.selection('input').set(strcat('as_roty_', int2str(obj.id)));
            rotx.set('pos', [0, 0, 0]);
            rotx.set('axis', [1, 0, 0]);
            rotx.set('rot',obj.theta_x);

            v = [obj.x, obj.y, obj.z];
            as = geom1.feature.create(strcat('as_', int2str(obj.id)),'Move');
            as.selection('input').set(strcat('as_rotx_', int2str(obj.id)));
            as.set('displ', v);
            as.set('createselection', 'on');

            %% Create geometry
            geom1.run();

            %% Create current terminal associated to the active site
            assel = mphgetselection(model.selection(strcat('geom1_', as.tag.char, '_dom')));
            asdom = assel.entities;
            model.physics('ec').create(strcat('term_', int2str(obj.id)), 'DomainTerminal', 3);
            model.physics('ec').feature(strcat('term_', int2str(obj.id))).label(strcat('term_', int2str(obj.id)));
            model.physics('ec').feature(strcat('term_', int2str(obj.id))).selection.set(asdom);

        end

        function model = assign_materials(obj, model)
            elecdom = mphgetselection(model.material('elec').selection()).entities;
            asdom = mphgetselection(model.material('tungsten').selection()).entities;

            elecsel = mphgetselection(model.selection(strcat('geom1_elec_', int2str(obj.id), '_dom')));
            elecdom = [elecdom, elecsel.entities];
            
            assel = mphgetselection(model.selection(strcat('geom1_as_', int2str(obj.id), '_dom')));
            asdom = [asdom, assel.entities];

            model.material('elec').selection.set(elecdom);
            model.material('tungsten').selection.set(asdom);
        end
        
        function model = set_active_site_current(obj, model, site_id, curr)
            if site_id ~= 1
                error("only 1 active site admitted for uNG")
            end
            model.physics('ec').feature(strcat('term_', int2str(obj.id))).set('I0', curr);
        end
        
        function [elec_sel, as_sel, encaps_sel] = get_selections(obj)
            elec_sel = {strcat('elec_', num2str(obj.id))};
            as_sel = {strcat('as_', int2str(obj.id))};
            encaps_sel = {};
        end

    end
end