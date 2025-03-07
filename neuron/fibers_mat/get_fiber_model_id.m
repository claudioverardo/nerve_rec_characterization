function fiber_model_id = get_fiber_model_id(fiber_model)
% Retrieve the internal id associated to a biophysical model of fibers.
% See build_fiber.py for the correspondent NEURON implementations.
%   
%   fiber_model_id = get_fiber_model_id(fiber_model)
% 
%   ---------------------------------------------------------------------------
%   fiber_model             id  myel    unmyel  ref (code origin)
%   ---------------------------------------------------------------------------
%   "ascent_sundt15"        1           x       Musselman et al 2021
%   "ascent_mrg"            6   x               Musselman et al 2021
%   ---------------------------------------------------------------------------
% 
% NOTE: also the opposite conversion id => string is supported.
    
    if ischar(fiber_model)
        fiber_model = string(fiber_model);
    end

    fiber_model_id = arrayfun(@(x) get_fiber_model_id_aux(x), fiber_model);

end

function fiber_model_id = get_fiber_model_id_aux(fiber_model)

    if isstring(fiber_model) || ischar(fiber_model)

        switch fiber_model
            case "ascent_sundt15"
                fiber_model_id = 1;
            case "ascent_mrg"
                fiber_model_id = 6;
            otherwise
                error("fiber_model not valid");
        end

    else

        switch fiber_model
            case 1
                fiber_model_id = "ascent_sundt15";
            case 6
                fiber_model_id = "ascent_mrg";
            otherwise
                error("fiber_model not valid");
        end

    end
end