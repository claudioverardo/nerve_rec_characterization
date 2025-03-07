classdef FascicleTopography
    % FascicleTopography Fascicle topography composed by fascicles with circular cross-section.
    % 
    % FascicleTopography Properties:
    %     length                        - length of the fascicles
    %     fascicles                     - geometry of the fascicles
    %     perineurium_type              - rule used to compute the thickness of the perineurium
    % 
    % FascicleTopography Properties (Dependent):
    %     n_fasc                        - number of fascicles
    % 
    % FascicleTopography Methods:
    %     FascicleTopography    - constructor
    %     add_to_model                  -
    %     add_fascicles_to_model        -
    %     add_ci_fascicles_to_model     -
    %     assign_materials              -
    %     get_perineurium_thickness     -
    % 
    % FascicleTopography Methods (Static):
    %     
    %     perineurium_thickness         - 

    properties
        % length - length of the fascicles (double) [m]
        length
        
        % fascicles - geometry of the fascicles (n_fasc x 3 double)
        % the i-th fascicle is described as follows:
        %     fascicles(i,:) = [xc, yc, r] 
        %     xc, yc: coordinates of the center of the fascicle cross-section [m]
        %     r: radius of the fascicle cross-section [m]
        fascicles

        % perineurium_type - rule used to compute the thickness of the perineurium 
        % available choices:
        %     "grinberg08" (default)
        %     "pelot20_human"
        %     "pelot20_pig"
        %     "pelot20_rat"
        %     (n_fasc x 1 double): the i-th entry defines the perineurium thickness [m] of the i-th fascicle
        %     fn handle: the input is the fascicle diameter [m] and the output the perineurium thickness [m]
        perineurium_type
    end

    properties (Dependent)
        % n_fasc - number of fascicles (int)
        n_fasc
    end
    
    methods

        function NFasc = get.n_fasc(obj)
            NFasc = size(obj.fascicles, 1);
        end

        function obj = FascicleTopography(params)
            % Instantiate object properties from struct of parameters
            % 
            % obj = FascicleTopography(params)

            obj.fascicles = params.fascicles;
            obj.length = params.nerve_extrusion_length;

            if isfield(params, 'perineurium_type')
                if isa(params.perineurium_type, 'numeric')
                    if length(params.perineurium_type) ~= obj.n_fasc
                        error("The perineurium thickness must be specified for each fascicle")
                    end
                end
                obj.perineurium_type = params.perineurium_type;
            else
                obj.perineurium_type = "grinberg08";
            end
        end
        
        function model = add_to_model(obj, model)
            % add fascicles to the geometry and perineurium contact impedance boundary conditions
            % 
            % model = add_to_model(model)

            model = add_fascicles_to_model(obj, model);
            model = add_ci_fascicles_to_model(obj, model);
        end
        
        function model = add_fascicles_to_model(obj, model)
            % model = add_fascicles_to_model(model)

            geom1 = model.geom('geom1');

            % create endoneurium
            endonames = {};
            for i =1:obj.n_fasc
                % Create an endoneurial compartment
                x = obj.fascicles(i, 1);
                y = obj.fascicles(i, 2);
                r = obj.fascicles(i, 3);
                
                endoaux(i) = geom1.feature.create(strcat('endoaux', int2str(i)), 'Cylinder');
                endoaux(i).set('h', obj.length);
                endoaux(i).set('x', x);
                endoaux(i).set('y', y);
                endoaux(i).set('r', r);
                endoaux(i).set('z', -obj.length/2);  % centers in z = 0
                
                endonames{end+1} = strcat('endoaux', int2str(i));
                
                % geom1.run;
            end
            
            % geom1.run;
            
            endogeom = geom1.feature.create('endogeom', 'Union');
            endogeom.selection('input').set(endonames);
            endogeom.set('createselection','on');
            
            endocopy = geom1.feature.create('endocopy', 'Copy');
            endocopy.selection('input').set('endogeom');
            
            geom1.run;

            epigeom = geom1.feature.create('epigeom','Difference');
            epigeom.selection('input').set('epifull');
            epigeom.selection('input2').set('endocopy');
            epigeom.set('createselection','on');
            
            geom1.run;
            
        end

        function model = add_ci_fascicles_to_model(obj, model)
            % model = add_ci_fascicles_to_model(model)
            
            geom1 = model.geom('geom1');
            N = 100;
            theta = linspace(0, 2*pi, N+1);
            % make perineural contact impedence
            for i = 1:obj.n_fasc
                % select the boundaries of each fascicle
                % NB. All the boundaries are selected (also top/bottom),
                % because it's easier and if the extrusion length is high
                % enough, it does not affect the solution.
                f(:, 1) = obj.fascicles(i, 1) + cos(theta) * obj.fascicles(i, 3);
                f(:, 2) = obj.fascicles(i, 2) + sin(theta) * obj.fascicles(i, 3);
                xmin = min(f(:, 1));
                ymin = min(f(:, 2));
                xmax = max(f(:, 1));
                ymax = max(f(:, 2));
                zmin = -obj.length/2;
                zmax = obj.length/2;
                delta = 0.01 * 1e-3;
                p0 = [xmin, ymin, zmin] - delta;
                p1 = [xmax, ymax, zmax] + delta;
                idx = mphselectbox(model, 'geom1', [p0', p1'], 'boundary');
                
                % deselect the top and bottom faces...
                % bottom
                p0 = [xmin, ymin, zmin] - delta;
                p1 = [xmax, ymax, zmin] + delta;
                idx_bottom = mphselectbox(model, 'geom1', [p0', p1'], 'boundary');
                
                % top
                p0 = [xmin, ymin, zmax] - delta;
                p1 = [xmax, ymax, zmax] + delta;
                idx_top = mphselectbox(model, 'geom1', [p0', p1'], 'boundary');
                
                idx = setdiff(idx, [idx_bottom, idx_top]);
                geom1.run;
                
                % calculate fascicle effective diameter
                % A = polyarea(f(:, 1), f(:, 2));
                % d = sqrt(A/pi)*2;
                % ds = 0.03 * d;
                ds = obj.get_perineurium_thickness(i);
                
                % create and apply contact impedance boundary conditions
                ci{i} = model.physics('ec').create(strcat('ci', int2str(i)), 'ContactImpedance', 2);
                ci{i}.selection.set(idx);
                ci{i}.set('ds', ds);
                ci{i}.set('sigmabnd_mat', 'userdef');
                ci{i}.set('epsilonrbnd_mat', 'userdef');
                ci{i}.set('sigmabnd', '0.0009[S/m]');
                ci{i}.set('epsilonrbnd', '80');
            end
            
            geom1.run;
            
        end
        
        function model = assign_materials(obj, model)
            % model = assign_materials(model)
            
            endosel = mphgetselection(model.selection('geom1_endogeom_dom'));
            endodom = endosel.entities;
            model.material('endo').selection.set(endodom);
            
        end

        function T_peri = get_perineurium_thickness(obj, idcs_fascicles)
            % Return perineurium thickness of fascicles [m].
            % 
            % NOTE: wrapper of FascicleTopography.perineurium_thickness().
            % 
            % T_peri = get_perineurium_thickness(idcs_fascicles)
            % 
            % PARAMETERS
            % ----------
            % idcs_fascicles (N x 1 double)
            %     ids of the fascicles
            % 
            % RETURNS
            % -------
            % T_peri (N x 1 double)
            %     thickness(es) of the perineurium of the fascicle(s) [m]
            
            if nargin == 1
                idcs_fascicles = 1:obj.n_fasc;
            end
            T_peri = zeros(length(idcs_fascicles),1);
            for i=1:length(idcs_fascicles)
                d = 2*obj.fascicles(idcs_fascicles(i),3);
                if isa(obj.perineurium_type, 'numeric')
                    T_peri(i) = obj.perineurium_type(idcs_fascicles(i));
                elseif isa(obj.perineurium_type, 'function_handle')
                    T_peri(i) = obj.perineurium_type(d);
                else
                    T_peri(i) = FascicleTopography.perineurium_thickness(obj.perineurium_type, d);
                end
            end
        end

    end

    methods (Static)

        function T_peri = perineurium_thickness(perineurium_type, d_fasc)
            % Return perineurium thickness of fascicles according to literature fits [m].
            % 
            % T_peri = perineurium_thickness(perineurium_type, d_fasc)
            % 
            % PARAMETERS
            % ----------
            % perineurium_type (string)
            %     "grinberg08", "pelot20_human", "pelot20_pig", "pelot20_rat" 
            % d_fasc (N x 1 double)
            %     (effective) diameter(s) of the fascicle(s) [m]
            % 
            % RETURNS
            % -------
            % T_peri (N x 1 double)
            %     thickness(es) of the perineurium of the fascicle(s) [m]

            if strcmp(perineurium_type, 'grinberg08')
                T_peri = 0.03 * d_fasc;
            elseif strcmp(perineurium_type, 'pelot20_human')
                T_peri = 0.03702 * d_fasc + 10.50e-6;
            elseif strcmp(perineurium_type, 'pelot20_pig')
                T_peri = 0.02547 * d_fasc + 3.440e-6;
            elseif strcmp(perineurium_type, 'pelot20_rat')
                T_peri = 0.01292 * d_fasc + 1.367e-6;
            else
                error("perineurium_type not valid")
            end
        end

    end
end