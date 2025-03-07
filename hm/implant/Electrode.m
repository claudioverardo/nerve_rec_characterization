classdef (Abstract) Electrode
    % Electrode Abstract class representing an electrode.
    
    properties
        % id - identifier of the electrode
        id
    end
    
    methods (Abstract)
        model = add_to_model(obj, model)
        model = assign_materials(obj, model)
        model = set_active_site_current(obj, model, site_id, curr)
        [elec_sel, as_sel, gnd_sel] = get_selections(obj)
    end
end