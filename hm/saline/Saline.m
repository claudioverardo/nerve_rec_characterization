classdef Saline
    % Saline Subclass of OuterNerve representing the surrounding saline bath.
    % 
    % Saline Properties:
    %     radius_ext                -
    %     topbot                    -
    %     
    % Saline Methods:
    %     Saline                    - constructor
    %     add_to_model              -
    %     assign_materials          -

    properties
        % radius_ext - radius of the cylindric domain of the saline bath (double) [m]
        radius_ext

        % topbot - longitudinal offset between the ends of the saline and the nerve domains (double) [m]
        % NOTE: set to zero by default
        topbot
    end
    
    methods

        function obj = Saline(params)
            obj.radius_ext = params.radius_ext;
            if ~isfield(params, 'topbot')
                params.topbot = 0;
            end
            obj.topbot = params.topbot;
        end

        function model = add_to_model(obj, model, nerve)
            geom1 = model.geom('geom1');
            
            disp('Adding saline bath')
            
            % create saline bath
            salfull = geom1.feature.create('salfull','Cylinder');
            salfull.set('h', nerve.length + obj.topbot);
            salfull.set('r', obj.radius_ext);
            salfull.set('pos', [0,0,-(nerve.length + obj.topbot)/2]);
            salfull.set('createselection', 'on');

            geom1.run()

            % set the ground boundary
            salsel = mphgetselection(model.selection(strcat('geom1_salfull_bnd')));
            salbnd = salsel.entities;
            model.physics('ec').feature.create('gnd', 'Ground', 2);
            model.physics('ec').feature('gnd').selection.set(salbnd);
        end
        
        function model = assign_materials(obj, model)
            salsel = mphgetselection(model.selection('geom1_salgeom_dom'));
            saldom = salsel.entities;
            model.material('sal').selection.set(saldom);
        end
    end
end