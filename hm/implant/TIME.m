classdef TIME < Electrode
    % TIME Subclass of Electrode representing a TIME electrode
    % (TIME: transverse intrafascicular multichannel electrode).
    % 
    % TIME Properties:
    %     id                        - (inherited from Electrode)
    % 
    %     ----------- electrode shaft parameters ----------
    %     l_shaft                   -
    %     w_shaft                   -
    %     h_shaft                   -
    % 
    %     ------------ active sites parameters ------------
    %     l_cc                      -
    %     n_as                      -
    %     d_as                      -
    %     h_as                      -
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
    %     ----------- encapsulation parameters ------------
    %     thick_encaps              -
    % 
    % TIME Properties (Dependent):
    %     n_sites                   -
    %     volume_as                 -
    %     
    % TIME Methods:
    %     TIME                      - constructor
    %     add_to_model              -
    %     assign_materials          -
    %     set_active_site_current   -
    %     get_selections            -

    properties
        % l_shaft - length of the electrode shaft (double) [m]
        l_shaft

        % w_shaft - width of the electrode shaft (double) [m]
        w_shaft

        % h_shaft - height/depth of the electrode shaft (double) [m]
        h_shaft
        
        % l_cc - distance between same face active sites (double) [m]
        l_cc
        
        % n_as - number of active sites per face (int)
        n_as

        % d_as - diameter of active sites (double) [m]
        d_as

        % h_as - depth of active sites (double) [m]
        h_as
        
        % x - origin x-location of the electrode reference frame (double) [m]
        x

        % y - origin y-location of the electrode reference frame (double) [m]
        y

        % z - origin z-location of the electrode reference frame (double) [m]
        z
        
        % theta_x - orientation around x-axis of the electrode reference frame (double) [deg]
        theta_x
        
        % theta_y - orientation around y-axis of the electrode reference frame (double) [deg]
        theta_y
        
        % theta_z - orientation around z-axis of the electrode reference frame (double) [deg]
        theta_z

        % thick_encaps - thickness of the encapsulation tissue (double) [deg]
        % set to zero (ie, no encapsulation) by default
        thick_encaps
    end
    
    properties (Dependent)
        % n_sites - total number of active sites (both faces) (int)
        %     n_sites = 2*n_as 
        n_sites

        % volume_as - volume of each active site (double) [m^3]
        volume_as
    end
    
    methods

        function NSites = get.n_sites(obj)
            NSites = 2 * obj.n_as;
        end

        function VolumeAs = get.volume_as(obj)
            VolumeAs = pi * (obj.d_as/2)^2 * obj.h_as;
        end
        
        function obj = TIME(params)
            obj.id = params.id;
            obj.l_shaft = params.l_shaft;
            obj.w_shaft = params.w_shaft;
            obj.h_shaft = params.h_shaft;
            obj.l_cc = params.l_cc;
            obj.n_as = params.n_as;
            obj.d_as = params.d_as;
            obj.h_as = params.h_as;

            obj.x = params.x;
            obj.y = params.y;
            obj.z = params.z;
            obj.theta_x = params.theta_x;
            obj.theta_y = params.theta_y;
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
            disp('Generating electrode geometry in electrode reference frame')
            %% Generate electrode body in electrode reference frame
            geom1 = model.geom('geom1');
            elec_full = geom1.create(strcat('elec_full_', int2str(obj.id)), 'Block');
            elec_full.set('size', [obj.l_shaft, obj.h_shaft, obj.w_shaft]);
            elec_full.set('base', 'center');
            elec_full.set('createselection', 'on');
            
            % geom1.run()
            
            %% Generate active sites in electrode reference frame
            as_names = {};
            for j = 1:obj.n_as
                as(j) = geom1.create(strcat('as_base_', int2str(obj.id), '_', int2str(j)), 'Cylinder');
                as(j).set('pos', [(j-1) * obj.l_cc - (obj.n_as-1) * obj.l_cc/ 2, -obj.h_shaft/2, 0]);
                as(j).set('r', obj.d_as/2);
                as(j).set('h', obj.h_as);
                as(j).set('axistype', 'y');
                as(j).set('createselection', 'on');
                
                % geom1.run();
                               
                as(obj.n_as + j) = geom1.create(strcat('as_base_', int2str(obj.id), '_', int2str(j + obj.n_as)), 'Cylinder');
                as(obj.n_as + j).set('r', obj.d_as/2);
                as(obj.n_as + j).set('h', obj.h_as);
                as(obj.n_as + j).set('axistype', 'y');
                as(obj.n_as + j).set('pos', [(j-1) * obj.l_cc - (obj.n_as-1) * obj.l_cc/ 2, obj.h_shaft/2 - obj.h_as, 0]);
                as(obj.n_as + j).set('createselection', 'on');

                % disp("(x,y) location of TIME's active sites")
                % disp((j-1) * obj.l_cc - (obj.n_as-1) * obj.l_cc/ 2)
                % disp(obj.h_shaft/2 - obj.h_as)
                
                % geom1.run();

                as_names = [as_names, strcat('as_base_', int2str(obj.id), '_', int2str(j)), strcat('as_base_', int2str(obj.id), '_', int2str(j + obj.n_as))];
            end
            
            elec = geom1.create(strcat('elec_base_', int2str(obj.id)), 'Difference');
            elec.selection('input').set(strcat('elec_full_', int2str(obj.id)));
            elec.selection('input2').set(as_names);
            elec.set("keepsubtract", true);
            
            % geom1.run()

            % Previous location of boundaries definition
            
            disp('Moving electrode geometry to laboratory reference frame')
            %% Apply roto-translation to electrode body
            % the resulting electrode body is tagged 'elec_[obj.id]'
            rotz = geom1.feature.create(strcat('elec_rotz_', int2str(obj.id)),'Rotate');
            rotz.selection('input').set(strcat('elec_base_', int2str(obj.id)));
            rotz.set('pos', [0, 0, 0]);
            rotz.set('axis', [0, 0, 1]);
            rotz.set('rot', obj.theta_z);
            % geom1.run();
            
            roty = geom1.feature.create(strcat('elec_roty_', int2str(obj.id)),'Rotate');
            roty.selection('input').set(strcat('elec_rotz_', int2str(obj.id)));
            roty.set('pos', [0, 0, 0]);
            roty.set('axis', [0, 1, 0]);
            roty.set('rot', obj.theta_y);
            % geom1.run();
            
            rotx = geom1.feature.create(strcat('elec_rotx_', int2str(obj.id)),'Rotate');
            rotx.selection('input').set(strcat('elec_roty_', int2str(obj.id)));
            rotx.set('pos', [0, 0, 0]);
            rotx.set('axis', [1, 0, 0]);
            rotx.set('rot',obj.theta_x);
            % geom1.run();
            
            v = [obj.x, obj.y, obj.z];
            electrode = geom1.feature.create(strcat('elec_', int2str(obj.id)),'Move');
            electrode.selection('input').set(strcat('elec_rotx_', int2str(obj.id)));
            electrode.set('displ', v);
            electrode.set('createselection', 'on');
            % geom1.run();
            
            %% Apply roto-translation to active sites
            % the resulting active sites is are tagged 'as_[obj.id]_[j]'
            rotz = geom1.feature.create(strcat('as_rotz_', int2str(obj.id)),'Rotate');
            rotz.selection('input').set(as_names);
            rotz.set('pos', [0, 0, 0]);
            rotz.set('axis', [0, 0, 1]);
            rotz.set('rot', obj.theta_z);
            % geom1.run();

            roty = geom1.feature.create(strcat('as_roty_', int2str(obj.id)),'Rotate');
            roty.selection('input').set(strcat('as_rotz_', int2str(obj.id)));
            roty.set('pos', [0, 0, 0]);
            roty.set('axis', [0, 1, 0]);
            roty.set('rot', obj.theta_y);
            % geom1.run();

            rotx = geom1.feature.create(strcat('as_rotx_', int2str(obj.id)),'Rotate');
            rotx.selection('input').set(strcat('as_roty_', int2str(obj.id)));
            rotx.set('pos', [0, 0, 0]);
            rotx.set('axis', [1, 0, 0]);
            rotx.set('rot',obj.theta_x);
            % geom1.run();

            v = [obj.x, obj.y, obj.z];
            as = geom1.feature.create(strcat('as_', int2str(obj.id)),'Move');
            as.selection('input').set(strcat('as_rotx_', int2str(obj.id)));
            as.set('displ', v);
            as.set('createselection', 'on');
            % geom1.run();

            %% Create encapsulation tissue
            if obj.thick_encaps > 0
                encaps_full = geom1.create(strcat('encaps_full_', int2str(obj.id)), 'Block');
                encaps_full.set('size', [obj.l_shaft + 2*obj.thick_encaps, obj.h_shaft + 2*obj.thick_encaps, obj.w_shaft + 2*obj.thick_encaps]);
                encaps_full.set('base', 'center');
                % geom1.run();
            
                v = [obj.x, obj.y, obj.z];
                encaps_moved = geom1.feature.create(strcat('encaps_moved_', int2str(obj.id)), 'Move');
                encaps_moved.selection('input').set(strcat('encaps_full_', int2str(obj.id)));
                encaps_moved.set('displ', v);
                % geom1.run();
    
                encaps = geom1.create(strcat('encaps_', int2str(obj.id)), 'Difference');
                encaps.selection('input').set(encaps_moved.tag.char);
                encaps.selection('input2').set({electrode.tag.char, as.tag.char});
                encaps.set("keepsubtract", true);
                encaps.set('createselection', 'on');
            end

            geom1.run();

            for j = 1:2*obj.n_as
                assel = mphgetselection(model.selection(strcat('geom1_as_base_', int2str(obj.id),'_', int2str(j), '_dom')));
                asdom = assel.entities;
                % model.physics('ec').create(strcat('cs_', int2str(obj.id),'_', int2str(j)), 'CurrentSource', 3);
                model.physics('ec').create(strcat('term_', int2str(obj.id),'_', int2str(j)), 'DomainTerminal', 3);
                model.physics('ec').feature(strcat('term_', int2str(obj.id),'_', int2str(j))).label(strcat('term_', int2str(obj.id),'_', int2str(j)));
                model.physics('ec').feature(strcat('term_', int2str(obj.id),'_', int2str(j))).selection.set(asdom);
            end

        end

        function model = assign_materials(obj, model)
            elecdom = mphgetselection(model.material('elec').selection()).entities;
            asdom = mphgetselection(model.material('as').selection()).entities;

            elecsel = mphgetselection(model.selection(strcat('geom1_elec_', int2str(obj.id), '_dom')));
            elecdom = [elecdom, elecsel.entities];
            
            assel = mphgetselection(model.selection(strcat('geom1_as_', int2str(obj.id), '_dom')));
            asdom = [asdom, assel.entities];

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
            % model.physics('ec').feature(strcat('cs_', int2str(obj.id),'_', int2str(site_id))).set('Qj', curr / obj.volume_as);
            model.physics('ec').feature(strcat('term_', int2str(obj.id),'_', int2str(site_id))).set('I0', curr);
        end
        
        function [elec_sel, as_sel, encaps_sel] = get_selections(obj)
            elec_sel = {strcat('elec_', num2str(obj.id))};
            as_sel = {};
            for k = 1:obj.n_sites
                as_sel{end+1} = strcat('as_base_', int2str(obj.id), '_', int2str(k));
            end
            if obj.thick_encaps>0
                encaps_sel = {strcat('encaps_', num2str(obj.id))};
            else
                encaps_sel = {};
            end
        end
        
    end
end