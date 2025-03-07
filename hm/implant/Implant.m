classdef Implant
    % Implant Implanted device composed by one or more electrodes.
    % 
    % WARNING: currently, only a single electrode per implant is fully supported 
    %          by the class HMModel (see note at the end of this help).
    % 
    % Implant Properties:
    %     implant_type              - type of the implant
    %     electrodes                - definition of the electrodes
    % 
    % Implant Properties (Dependent):
    %     n_elecs                   - number of electrodes
    %
    % Implant Methods:
    %     Implant                   - constructor
    %     add_to_model              -
    %     assign_materials          -
    %     set_active_site_current   -

    properties
        % implant_type - type of the implant (string)
        % available choices: 'Cortec_Cuff', 'TIME', 'uNG'
        implant_type

        % electrodes - definition of the electrodes (n_elecs x 1 Electrode subclass)
        electrodes
    end
    
    properties (Dependent)
        n_elecs
    end
    
    methods

        function NElecs = get.n_elecs(obj)
            NElecs = length(obj.electrodes);
        end

        function obj = Implant(params)

            obj.implant_type = params.implant_type;
                
            for j = 1:params.n_elecs
                params.electrodes(j).id = params.id_elecs(j);

                switch obj.implant_type
                    case 'Cortec_Cuff'
                        electrodes(j) = CortecCuff(params.electrodes(j));
                    case 'TIME'
                        electrodes(j) = TIME(params.electrodes(j));
                    case 'uNG'
                        electrodes(j) = uNG(params.electrodes(j));
                    otherwise
                        error("Implant '"+obj.implant_type+"' not valid");
                end
    
            end
    
            obj.electrodes = electrodes;
        end
        
        function model = add_to_model(obj, model)
            for i = 1:obj.n_elecs
                disp(strcat("# Adding electrode n. ", int2str(i), " of implant ", obj.implant_type))
                model = obj.electrodes(i).add_to_model(model);
            end
        end
        
        function model = assign_materials(obj, model)
            for i = 1:obj.n_elecs
                model = obj.electrodes(i).assign_materials(model);
            end
        end
        
        function model = set_active_site_current(obj, elec_id, site_id, curr)
            model = obj.electrodes(elec_id).set_active_site_current(site_id, curr);
        end
    end

end