classdef FiberTopography < handle
    % FiberTopography Population of 2D-extruded fibers along with their groupings.
    % 
    % FiberTopography Properties:
    %     topo_func                     - definition of the fibers (position, diameter, biophysics)
    %     dir_topo_func                 - directory where topo_func.mat is stored
    % 
    % FiberTopography Properties (Dependent):
    %     n_fibers                      - number of fibers
    %
    % FiberTopography Methods:
    %     FiberTopography               - constructor
    %     get_node_locs                 -
    %     get_SUAP_template             -
    %     get_slice_topo_func           -
    
    properties
        % topo_func - definition of the fibers (n_fibers x 6 double)
        % the i-th row defines the property of the i-th fiber as follows:
        %   topo_func(i,1) = id                 fiber identifier
        %   topo_func(i,2) = x                  fiber x-position in the 2D cross-section of the nerve [m]
        %   topo_func(i,3) = y                  fiber y-position in the 2D cross-section of the nerve [m]
        %   topo_func(i,4) = fiber_model_id     fiber biophysical model identifier *
        %   topo_func(i,5) = d_fiber            fiber diameter [um]
        %   topo_func(i,6) = z_shift            fiber longitudinal shift [-] **
        % 
        % * see get_fiber_model_id() for the available models
        % ** expressed as a fraction of:
        %    - internodal length for myelinated fibers
        %    - discretization length for unmyelinated fibers
        topo_func

        % dir_topo_func - directory where to store topo_func.mat (string)
        % NOTE 1: to avoid storing the topo_func.mat, set it to []
        % NOTE 2: if part of an HMModel object, it is set equal to HMModel.dir_model 
        dir_topo_func
    end

    properties (Dependent)
        % n_fibers - number of fibers (int)
        n_fibers 
    end
    
    methods

        function n_fibers = get.n_fibers(obj)
            n_fibers = size(obj.topo_func, 1);
        end

        function set.topo_func(obj, topo_func) 
            obj.topo_func = topo_func;

            % The set method is invoked when the class is loaded from a mat file to initialize
            % all the class properties to [] and then to copy their values from the mat file.
            % In such cases, dir_topo_func is empty and therefore a copy of topo_func is created 
            % in the current directory. The if-statement circumvent such situations.
            if ~isempty(obj.dir_topo_func)
                save(fullfile(obj.dir_topo_func,"topo_func.mat"), 'topo_func');
            end
        end

        function obj = FiberTopography(params)
            % Instantiate object properties from struct of parameters
            % 
            % obj = FiberTopography(params)

            obj.dir_topo_func = params.dir_model;
            obj.topo_func = params.topo_func;
        end

        function [node_locs, node_locs_metadata] = get_node_locs(obj, id_fibers, nerve_length, varargin)
            % Generate 3D locations of fiber nodes (compartments) from fibers in the 2D functional topography.
            % An adaptive method is used to find the number of nodes corresponding to the largest fiber length
            % that is smaller than the nerve length.
            % 
            % [node_locs, node_locs_metadata] = get_node_locs(id_fibers, nerve_length, varargin)
            % 
            % PARAMETERS
            % ----------
            % id_fibers (n_fibers x 1 int)
            %     ids of the fibers 
            % nerve_length (double)
            %     length of the nerve
            % 
            % NAME, VALUE ARGS
            % ----------------
            % node_length_unmyel (double)
            %     spatial discretization step for unmyelinated fibers (50/6 um by default)
            % log_step (int)
            %     control amount of log information in the command window
            % n_internodes_myel_step (int)
            %     step of # internodes used with myelinated fibers to find the # nodes to fit the nerve length 
            %     IMPORTANT! EVEN NUMBER TO HAVE A RANVIER NODE @ z=0 (before z-shift)
            % n_nodes_unmyel_step (int)
            %     step of # nodes used with unmyelinated fibers to find the # nodes to fit the nerve length 
            %     IMPORTANT! ODD NUMBER TO HAVE A NODE @ z=0 (before z-shift)
            %
            % RETURNS
            % -------
            % node_locs (n_fibers x 3 array)
            %     content of locs.txt files (see HMModel.export_node_locs)
            % node_locs_metadata (struct)
            %     content of locs_meta.mat files (see HMModel.export_node_locs)

            if ~exist("id_fibers",'var') || isempty(id_fibers)
                id_fibers = obj.topo_func(:,1);
            end
            
            % input parser (hyperparameters used to generate fibers)
            p = inputParser;
            addParameter(p, 'n_internodes_myel_step', 100);
            addParameter(p, 'n_nodes_unmyel_step', 1251);
            addParameter(p, 'node_length_unmyel', 50/6);    % 8.3333 um - default in ASCENT
            addParameter(p, 'log_step', 50);
            p.KeepUnmatched = true;
            parse(p, varargin{:});
            
            % parse function parameters
            n_internodes_myel_step = p.Results.n_internodes_myel_step;
            n_nodes_unmyel_step = p.Results.n_nodes_unmyel_step;
            node_length_unmyel = p.Results.node_length_unmyel;
            log_step = p.Results.log_step;
            % args_unmatched = p.Unmatched;
        
            % select the properties of the input fibers
            topo_func_sel = obj.get_slice_topo_func(id_fibers);
            n_fibers_sel = size(topo_func_sel,1);

            % [x,y,z] 3D positions of fiber nodes
            node_locs = [];
            
            % metadata of the generation of fiber nodes
            z_shift_m = zeros(n_fibers_sel,1);
            n_nodes = zeros(n_fibers_sel,1);
            n_internodes_myel = zeros(n_fibers_sel,1);
            n_nodes_unmyel = zeros(n_fibers_sel,1);
            locs_id_fiber = [];
            locs_mask_ranvier = false(0);
            locs_idcs_inside_nerve = [];
        
            % sweep over selected fibers
            disp("i/tot) id x y model_id diam z_shift");
            changed_fiber_mod = true;
            for i=1:n_fibers_sel
        
                id_fiber = topo_func_sel(i,1);
                x = topo_func_sel(i,2);
                y = topo_func_sel(i,3);
                if i>1
                    changed_fiber_mod = fiber_model_id ~= topo_func_sel(i,4) || d_fiber ~= topo_func_sel(i,5);
                end
                fiber_model_id = topo_func_sel(i,4);
                fiber_model = get_fiber_model_id(fiber_model_id);
                d_fiber = topo_func_sel(i,5);
                z_shift = topo_func_sel(i,6);
                
                if i==1 || mod(i,log_step)==0 || i==n_fibers_sel
                    disp(i+"/"+n_fibers_sel+") "+id_fiber+" "+num2str(x)+" "+num2str(y)+" "+fiber_model+" "+d_fiber+" "+z_shift);
                end
    
                % =========================== fiber node z-positions ===========================
                % before z-shift and w/o resize of fiber length to nerve length

                if changed_fiber_mod

                    % ========== from matlab surrogate of build_fiber.py
                    % it checks that the generated virtual fiber is longer than the fem nerve (for any z-shift)
                    % "virtual" fiber since only the geometry is built, not the biophysics
                    n_internodes_myel_current = 0;
                    n_nodes_unmyel_current = 0;
                    fiber_params = [];
                    if check_fiber_myelinated(fiber_model)
                        fiber_params.fiber_model = fiber_model;
                        fiber_params.diameter = d_fiber;
                        fiber_params.n_internodes = 0; % keep even
    
                        fiber_longer_than_nerve = false;
                        while ~fiber_longer_than_nerve
                            fiber_params.n_internodes = fiber_params.n_internodes + n_internodes_myel_step;  % keep even
                            [fiber_pos_z, fiber_mask_ranvier, fiber_morpho] = get_fiber_pos_z(fiber_params); % mm
                            fiber_characteristic_length = fiber_morpho.length_internode; % um
                            fiber_longer_than_nerve = (fiber_pos_z(end)-fiber_pos_z(1))*1e-3 - 0.5*fiber_characteristic_length*1e-6 > nerve_length;
                        end
                        n_internodes_myel_current = fiber_params.n_internodes;
                    
                    elseif check_fiber_unmyelinated(fiber_model)
                        fiber_params.fiber_model = fiber_model;
                        fiber_params.node_length = node_length_unmyel;
                        fiber_params.n_nodes = 1; % keep odd
    
                        fiber_longer_than_nerve = false;
                        while ~fiber_longer_than_nerve
                            fiber_params.n_nodes = fiber_params.n_nodes + n_nodes_unmyel_step - 1; % keep odd
                            [fiber_pos_z, fiber_mask_ranvier, ~] = get_fiber_pos_z(fiber_params);       % mm
                            fiber_characteristic_length = fiber_params.node_length; % um
                            fiber_longer_than_nerve = (fiber_pos_z(end)-fiber_pos_z(1))*1e-3 - 0.5*fiber_characteristic_length*1e-6 > nerve_length;
                        end
                        n_nodes_unmyel_current = fiber_params.n_nodes;

                    else
                        error("not valid fiber_id");
                    end
                end

                n_internodes_myel(i) = n_internodes_myel_current;
                n_nodes_unmyel(i) = n_nodes_unmyel_current;
            
                % ==================== fiber node z-positions [INSIDE NERVE] =====================
                % apply z-shift and keep only the fiber nodes with positions that are inside the fem nerve
                % NOTE: the z-shift is handled differently in myelinated/unmyelinated fibers (see docs)

                z_shift_m(i) = fiber_characteristic_length*1e-6*z_shift; % m
                node_locs_z_i = fiber_pos_z*1e-3 + z_shift_m(i);
                
                idcs_inside_nerve = find(...
                    (node_locs_z_i >= -nerve_length/2) & ...
                    (node_locs_z_i <= nerve_length/2));
                
                node_locs_z_i = node_locs_z_i(idcs_inside_nerve);
                node_locs_mask_ranvier_i = fiber_mask_ranvier(idcs_inside_nerve);
                n_nodes(i) = length(node_locs_z_i);
                
                % store results
                node_locs = [node_locs; repmat([x y],n_nodes(i),1) node_locs_z_i(:)];
                locs_id_fiber = [locs_id_fiber; repmat(id_fiber,n_nodes(i),1)];
                locs_mask_ranvier = [locs_mask_ranvier; node_locs_mask_ranvier_i(:)];
                locs_idcs_inside_nerve = [locs_idcs_inside_nerve; idcs_inside_nerve(:)];
        
                % NOTE: not so important to keep track of "locs_idcs_inside_nerve". It was used in
                % previous versions of the code to compute the SUAPs from excitation templates.
                % Now, another heuristic is used.
            end

            node_locs_metadata.locs_id_fiber = locs_id_fiber;
            node_locs_metadata.locs_mask_ranvier = locs_mask_ranvier;
            node_locs_metadata.n_nodes = n_nodes;
            node_locs_metadata.z_shift_m = z_shift_m;
            node_locs_metadata.node_length_unmyel = node_length_unmyel;
            node_locs_metadata.n_internodes_myel = n_internodes_myel;
            node_locs_metadata.n_nodes_unmyel = n_nodes_unmyel;
            node_locs_metadata.n_internodes_myel_step = n_internodes_myel_step; % DEPRECATED
            node_locs_metadata.n_nodes_unmyel_step = n_nodes_unmyel_step;       % DEPRECATED
            node_locs_metadata.locs_idcs_inside_nerve = locs_idcs_inside_nerve; % DEPRECATED
        
        end

        function SUAP_data = get_SUAP_template(obj, id_fibers, nerve_length, node_locs, z_shift_m, locs_id_fiber, locs_mask_ranvier, LFM, dir_fiber_templates, varargin)
            % Compute SUAPs of fibers starting from templates of fiber excitations.
            % 
            % SUAP_data = get_SUAP_template(id_fibers, nerve_length, node_locs, z_shift_m, locs_id_fiber, locs_mask_ranvier, LFM, dir_fiber_templates, varargin)
            %
            % PARAMETERS
            % ----------
            % id_fibers (n_fibers x 1)
            %     ids of the fibers 
            % nerve_length (double)
            %     length of the nerve
            % node_locs
            %     see get_node_locs()
            % z_shift_m
            %     see get_node_locs()
            % locs_id_fiber
            %     see get_node_locs()
            % locs_mask_ranvier
            %     see get_node_locs()
            % LFM (n_locs x 1 double)
            %     LFM of the active site. The LFM is used to compute the transresistance 
            %     of fiber nodes according to the reciprocity theorem of electromagnetism.
            %     NOTE: NOT (n_locs x 4) as in LFM_#elec_#as.txt files
            % dir_fiber_templates (string)
            %     path to the directory containing the excitation templates
            % 
            % NAME, VALUE ARGS
            % ----------------
            % v_thresh (double)
            %     -40 [mV] by default
            %     threshold for the membrane potential to detect the spike initiation
            % log_step (int)
            %     control amount of log information in the command window
            % 
            % RETURNS
            % -------
            % SUAP_data (struct)
            %     content of SUAP_#elec_#as.mat files (see HMModel.export_SUAP_template)

            % input parser
            p = inputParser;
            addParameter(p, 'v_thresh', -40);
            addParameter(p, 'log_step', 50);
            parse(p, varargin{:});
            
            % parse function parameters
            v_thresh = p.Results.v_thresh;
            log_step = p.Results.log_step;

            % legacy
            align_init_spike = false;

            % ================= retrieve the fiber topography of the nerve ================
            % select the properties of the input fibers
            [topo_func_sel, topo_func_sel_idcs] = obj.get_slice_topo_func(id_fibers);
            n_fibers_sel = size(topo_func_sel,1);

            % ======================== compute SUAPs of all fibers ========================
            disp("i/tot) id x y model_id diam z_shift template_type");
            SUAP = cell(n_fibers_sel,1);
            SUAP_t = cell(n_fibers_sel,1);
            SUAP_dt = zeros(n_fibers_sel,1);
            SUAP_idx_max = zeros(n_fibers_sel,1);
            SUAP_idx_min = zeros(n_fibers_sel,1);
            SUAP_template_type = strings(n_fibers_sel,1);
            SUAP_CV = zeros(n_fibers_sel,1);
            SUAP_align_init_spike = false(n_fibers_sel,1);
            SUAP_t0 = zeros(n_fibers_sel,1);
            SUAP_t0_idx = zeros(n_fibers_sel,1);
            changed_fiber_mod = true;
            for i=1:n_fibers_sel
        
                % ======================= retrieve the fiber properties =======================
                % needed to retrieve the correct excitation template / data from the LFM
                id_fiber = topo_func_sel(i,1);
                x = topo_func_sel(i,2);
                y = topo_func_sel(i,3);
                if i>1
                    changed_fiber_mod = fiber_model_id ~= topo_func_sel(i,4) || d_fiber ~= topo_func_sel(i,5);
                end
                fiber_model_id = topo_func_sel(i,4);
                fiber_model = get_fiber_model_id(fiber_model_id);
                d_fiber = topo_func_sel(i,5);
                z_shift = topo_func_sel(i,6);

                % ========== retrieve the LFM/locs values at the fiber compartments ==========
                locs_fiber_idcs = find(locs_id_fiber == id_fiber);
                locs_fiber_pos_z = node_locs(locs_fiber_idcs,3);
                LFM_fiber = LFM(locs_fiber_idcs,1);
                    
                % positions of Ranvier nodes
                locs_fiber_mask_ranvier = locs_mask_ranvier(locs_fiber_idcs);
                locs_fiber_idcs_ranvier = find(locs_fiber_mask_ranvier);
                locs_fiber_n_ranvier = length(locs_fiber_idcs_ranvier);
    
                % ================== retrieve the fiber excitation template ===================
                if changed_fiber_mod
                    fiber_template_basepath = fullfile(dir_fiber_templates,fiber_model+"/template_fiber_d_"+num2str(d_fiber,"%.2f"));
                    if ~isfile(fiber_template_basepath+".mat")
                        % legacy reasons
                        fiber_template_basepath = fullfile(dir_fiber_templates,fiber_model+"/template_fiber_d_"+num2str(d_fiber,"%.1f")); 
                    end
                    % template data
                    % NOTE: fiber_t_rec / fiber_v_rec used only if fiber_template_type = "internode-CV"
                    load(fiber_template_basepath+".mat", 'fiber_t_rec', 'fiber_v_rec', 'fiber_itot_rec');
                    % template metadata
                    % NOTE: fiber_CV used only if fiber_template_type = "internode-CV"
                    %       fiber_pos_z used only if fiber_template_type = "full-fiber"
                    load(fiber_template_basepath+"_meta.mat", 'fiber_template_type', 'fiber_area', 'fiber_pos_z', 'fiber_CV', 'fiber_dt'); 
                    fiber_pos_z = fiber_pos_z*1e-3; % m
                end
                
                % show log
                if i==1 || mod(i,log_step)==0 || i==n_fibers_sel
                    disp(i+"/"+n_fibers_sel+") "+id_fiber+" "+num2str(x)+" "+num2str(y)+" "+fiber_model+" "+d_fiber+" "+z_shift+" "+fiber_template_type);
                end
    
                % =================== CASE 1: template type "internode-CV" ====================
                % In the excitation templates are stored the variables (v, itot) of a single fiber internode.
                if strcmp(fiber_template_type, 'internode-CV')

                    % store the excitation template in separated variables
                    % consistently w/ the "full-fiber" case, the excitation template along all the nerve
                    % should be contained in fiber_t_rec, fiber_area, fiber_v_rec, fiber_itot_rec 
                    if changed_fiber_mod
                        fiber_t_rec_internode = fiber_t_rec;
                        fiber_v_rec_internode = fiber_v_rec;
                        fiber_itot_rec_internode = fiber_itot_rec;
                        fiber_area_internode = fiber_area;
                    end

                    % ========================= myelinated fiber
                    if locs_fiber_n_ranvier > 0

                        fiber_internode_n_nodes = size(fiber_itot_rec_internode,1);
                        fiber_internode_length = locs_fiber_pos_z(locs_fiber_idcs_ranvier(2))-locs_fiber_pos_z(locs_fiber_idcs_ranvier(1));

                        % data structures to store the reconstructed itot
                        fiber_t_rec_internode_duration = fiber_t_rec_internode(end)-fiber_t_rec_internode(1); % [ms]
                        fiber_t_rec_duration = (locs_fiber_n_ranvier+1)*fiber_internode_length / fiber_CV * 1e+3 + fiber_t_rec_internode_duration; % [ms]
                        fiber_t_rec_n_samples = ceil(fiber_t_rec_duration / fiber_dt);
                        fiber_t_rec = (0:fiber_t_rec_n_samples-1)*fiber_dt;
                        fiber_itot_rec = zeros((locs_fiber_n_ranvier+1)*fiber_internode_n_nodes, fiber_t_rec_n_samples);
                        % fiber_v_rec = zeros((locs_fiber_n_ranvier+1)*fiber_internode_n_nodes, fiber_t_rec_n_samples);
    
                        % faster than interp1()
                        F_interp = {};
                        for k=1:fiber_internode_n_nodes
                            F_interp{k} = griddedInterpolant(fiber_t_rec_internode, fiber_itot_rec_internode(k,:));
                        end
                    
                        % build itot(t) at each fiber nodes by time-shifting the excitation template according to
                        % the conduction velocity and the locations of the nodes
                        for n=1:locs_fiber_n_ranvier+1
                            
                            delay_n_internode = (n-1)*fiber_internode_length / fiber_CV * 1e+3;
                            delay_n_internode_n_samples = floor(delay_n_internode/fiber_dt) + 1;
                            t_n_internode_idcs = delay_n_internode_n_samples + (0:length(fiber_t_rec_internode)-1);
                            t_n_internode_idcs = t_n_internode_idcs(t_n_internode_idcs<=length(fiber_t_rec));
                            
                            % itot_n_internode = fiber_itot_rec_internode; % bad
                            itot_n_internode = zeros(fiber_internode_n_nodes,length(fiber_t_rec_internode));
                            % v_n_internode = zeros(fiber_internode_n_nodes,length(fiber_t_rec_internode));
                            for k=1:fiber_internode_n_nodes
                                % itot_n_internode(k,:) = interp1(fiber_t_rec_internode + delay_n_internode, fiber_itot_rec_internode(k,:), fiber_t_rec(t_n_internode_idcs));
                                itot_n_internode(k,:) = F_interp{k}(fiber_t_rec(t_n_internode_idcs)-delay_n_internode);
                                % v_n_internode(k,:) = interp1(fiber_t_rec_internode + delay_n_internode, fiber_v_rec_internode(k,:), fiber_t_rec(t_n_internode_idcs));
                            end
                            itot_n_internode(isnan(itot_n_internode)) = 0;
                            fiber_itot_rec(((n-1)*fiber_internode_n_nodes+1):n*fiber_internode_n_nodes,t_n_internode_idcs) = itot_n_internode;
                            % fiber_v_rec(((n-1)*fiber_internode_n_nodes+1):n*fiber_internode_n_nodes,t_n_internode_idcs) = v_n_internode;
                            % fiber_v_rec(((n-1)*fiber_internode_n_nodes+1):n*fiber_internode_n_nodes,fiber_t_rec<fiber_t_rec(t_n_internode_idcs(1))) = fiber_v_rec_internode(1);
                            % fiber_v_rec(((n-1)*fiber_internode_n_nodes+1):n*fiber_internode_n_nodes,fiber_t_rec>fiber_t_rec(t_n_internode_idcs(end))) = fiber_v_rec_internode(end);
                            
                        end

                        % remove nodes at start/end of nerve that are not inside the nerve
                        % (itot has been generated as non-fractional time-shifted replicas of the template)
                        n_nodes_neglect_start = fiber_internode_n_nodes - (locs_fiber_idcs_ranvier(1)-1);
                        n_nodes_neglect_end = fiber_internode_n_nodes - (length(locs_fiber_idcs)-locs_fiber_idcs_ranvier(end)+1);
                        fiber_itot_rec = fiber_itot_rec((n_nodes_neglect_start+1):(end-n_nodes_neglect_end),:);
                        % fiber_v_rec = fiber_v_rec((n_nodes_neglect_start+1):(end-n_nodes_neglect_end+1),:);

                        % time instant when v > v_thresh @ z = beginning of the nerve
                        % (detect+shifted from fiber_v_rec_internode to avoid generating also fiber_v_rec)
                        [~, fiber_t_rec_internode_t0_idx] = find(fiber_v_rec_internode(1,:) > v_thresh);
                        fiber_t_rec_internode_t0_idx = fiber_t_rec_internode_t0_idx(1);
                        fiber_t_rec_internode_t0 = fiber_t_rec_internode(fiber_t_rec_internode_t0_idx);
                        fiber_t_rec_t0 = fiber_t_rec_internode_t0 + (fiber_internode_length - (locs_fiber_pos_z(locs_fiber_idcs_ranvier(1))-(-nerve_length/2))) / fiber_CV * 1e+3;
                        fiber_t_rec_t0_idx = round(fiber_t_rec_t0/fiber_dt);                        
                    
                        % for consistency with the "full-fiber"
                        fiber_area = repmat(fiber_area_internode, 1, locs_fiber_n_ranvier-1);
                        if n_nodes_neglect_start < fiber_internode_n_nodes
                            fiber_area = [fiber_area_internode((n_nodes_neglect_start+1):end) fiber_area];
                        end
                        if n_nodes_neglect_end < fiber_internode_n_nodes
                            fiber_area = [fiber_area fiber_area_internode(1:(end-n_nodes_neglect_end))];
                        end

                    % ========================= unmyelinated fiber
                    else

                        % data structures to store the reconstructed itot
                        fiber_t_rec_internode_duration = fiber_t_rec_internode(end)-fiber_t_rec_internode(1); % [ms]
                        fiber_t_rec_duration = (locs_fiber_pos_z(end)-locs_fiber_pos_z(1)) / fiber_CV * 1e+3 + fiber_t_rec_internode_duration; % [ms]
                        fiber_t_rec_n_samples = ceil(fiber_t_rec_duration / fiber_dt);
                        fiber_t_rec = (0:fiber_t_rec_n_samples-1)*fiber_dt;
                        fiber_itot_rec = zeros(length(locs_fiber_idcs), fiber_t_rec_n_samples);
                        % fiber_v_rec = zeros(length(locs_fiber_idcs), fiber_t_rec_n_samples);
    
                        % faster than interp1()
                        F_interp = griddedInterpolant(fiber_t_rec_internode,fiber_itot_rec_internode);

                        % build itot(t) at each fiber nodes by time-shifting the excitation template according to
                        % the conduction velocity and the locations of the nodes
                        for n=1:length(locs_fiber_idcs)
                            delay_n_node = (locs_fiber_pos_z(n)-locs_fiber_pos_z(1))/ fiber_CV * 1e+3;
                            delay_n_node_n_samples = floor(delay_n_node/fiber_dt) + 1;
                            t_n_node_idcs = delay_n_node_n_samples + (0:length(fiber_t_rec_internode)-1); 
                            t_n_node_idcs = t_n_node_idcs(t_n_node_idcs<=length(fiber_t_rec));
                            % itot_n_node = fiber_itot_rec_internode; % bad
                            % itot_n_node = interp1(fiber_t_rec_internode + delay_n_node, fiber_itot_rec_internode, fiber_t_rec(t_n_node_idcs));
                            itot_n_node = F_interp(fiber_t_rec(t_n_node_idcs)-delay_n_node);
                            itot_n_node(isnan(itot_n_node)) = 0;
                            fiber_itot_rec(n,t_n_node_idcs) = itot_n_node;
                            % v_n_node = interp1(fiber_t_rec_internode + delay_n_node, fiber_v_rec_internode, fiber_t_rec(t_n_node_idcs));
                            % fiber_v_rec(n,t_n_node_idcs) = v_n_node;
                            % fiber_v_rec(n,fiber_t_rec<fiber_t_rec(t_n_node_idcs(1))) = fiber_v_rec_internode(1);
                            % fiber_v_rec(n,fiber_t_rec>fiber_t_rec(t_n_node_idcs(end))) = fiber_v_rec_internode(end);
                        end
                        
                        % time instant when v > v_thresh @ z = beginning of the nerve
                        % (detect+shifted from fiber_v_rec_internode to avoid generating also fiber_v_rec)
                        [~, fiber_t_rec_internode_t0_idx] = find(fiber_v_rec_internode(1,:) > v_thresh);
                        fiber_t_rec_internode_t0_idx = fiber_t_rec_internode_t0_idx(1);
                        fiber_t_rec_internode_t0 = fiber_t_rec_internode(fiber_t_rec_internode_t0_idx);
                        fiber_t_rec_t0 = fiber_t_rec_internode_t0 - (locs_fiber_pos_z(1)-(-nerve_length/2)) / fiber_CV * 1e+3;
                        fiber_t_rec_t0_idx = round(fiber_t_rec_t0/fiber_dt);

                        % for consistency with the "full-fiber" case
                        node_length_unmyel = (locs_fiber_pos_z(2)-locs_fiber_pos_z(1))*1e+6; % [um]
                        fiber_node_length_unmyel = fiber_area_internode / (pi * d_fiber); % [um]
                        if abs(node_length_unmyel - fiber_node_length_unmyel) > 1e-9
                            warning("Spatial discretization in FEM ("+node_length_unmyel+" um) and template ("+fiber_node_length_unmyel+" um) are different");
                        end
                        fiber_area = pi * d_fiber * node_length_unmyel;
                        fiber_area = repmat(fiber_area, 1, length(locs_fiber_idcs));
                        % fiber_area = repmat(fiber_area_internode, 1, length(locs_fiber_idcs));

                    end

                    % ============ DEBUGGING PLOTS - fiber_itot_rec reconstruction ============
                    % figure;
                    % imagesc('XData', fiber_t_rec, 'YData', locs_fiber_pos_z, 'CData', fiber_itot_rec);
                    % % imagesc('XData', fiber_t_rec, 'YData', locs_fiber_pos_z, 'CData', fiber_v_rec);
                    % colormap(gray);
                    % xlabel('Time [ms]');
                    % ylabel('Longitudinal position [mm]');
                    % h_cb = colorbar;
                    % h_cb.Label.String = 'itot';
                    % ==========================================================================

                    % for consistency with the "full-fiber" case
                    LFM_fiber_padding = LFM_fiber;

                % ===================== CASE 2: template type "full-fiber" ====================
                % In the excitation templates are stored the variables (v, itot) of all the fiber compartments.
                % The fibers in the excitation templates are assumed longer than the fem nerve and w/o z-shift.
                elseif strcmp(fiber_template_type, 'full-fiber')

                    % Apply z-shift to the positions of the excitation template
                    fiber_pos_z_shifted = fiber_pos_z + z_shift_m(topo_func_sel_idcs(i)); % m
        
                    % Add zero-padding to the LFM outside the fem nerve to account for
                    % the mismatch w/ the fiber length in the excitation template
                    idcs_nodes_neuron2fem = find( ...
                        (fiber_pos_z_shifted >= min(locs_fiber_pos_z)-1e-12) & ...
                        (fiber_pos_z_shifted <= max(locs_fiber_pos_z)+1e-12));
                    LFM_fiber_padding = zeros(length(fiber_pos_z_shifted),1);
                    LFM_fiber_padding(idcs_nodes_neuron2fem) = LFM_fiber;

                    % ========== DEBUGGING PLOTS - check LFM/template locs alignment ==========
                    % figure; 
                    % plot(fiber_pos_z_shifted, ones(size(fiber_pos_z_shifted)), 'o'); 
                    % hold on; 
                    % plot(locs_fiber_pos_z, 2*ones(size(locs_fiber_pos_z)), 'o');
                    % plot(fiber_pos_z_shifted(idcs_nodes_neuron2fem), 3*ones(size(idcs_nodes_neuron2fem)), 'o');
                    % ylim([0 4]);
                    % ==========================================================================

                    % time instant when v > v_thresh @ z = beginning of the nerve
                    % myelinated fiber
                    if locs_fiber_n_ranvier > 0
                        [~, fiber_t_rec_t0_idx] = find(fiber_v_rec(idcs_nodes_neuron2fem(locs_fiber_idcs_ranvier(1)),:) > v_thresh);
                        fiber_t_rec_t0_idx = fiber_t_rec_t0_idx(1);
                        fiber_t_rec_t0 = fiber_t_rec(fiber_t_rec_t0_idx);
                        fiber_t_rec_t0 = fiber_t_rec_t0 - (locs_fiber_pos_z(locs_fiber_idcs_ranvier(1))-(-nerve_length/2)) / fiber_CV * 1e+3;
                    % unmyelinated fiber
                    else
                        [~, fiber_t_rec_t0_idx] = find(fiber_v_rec(idcs_nodes_neuron2fem(1),:) > v_thresh);
                        fiber_t_rec_t0_idx = fiber_t_rec_t0_idx(1);
                        fiber_t_rec_t0 = fiber_t_rec(fiber_t_rec_t0_idx);
                        fiber_t_rec_t0 = fiber_t_rec_t0 - (locs_fiber_pos_z(1)-(-nerve_length/2)) / fiber_CV * 1e+3;

                        warning("TODO: handle the case length_node_unmyel in FEM =/= in template");
                    end
                    fiber_t_rec_t0_idx = round(fiber_t_rec_t0/fiber_dt);

                else
                    error('template_type not valid')

                end

                % set t=0 of the excitation template when v > v_thresh (if required)
                if align_init_spike
                    fiber_t_rec = fiber_t_rec(:,fiber_t_rec_t0_idx:end) - fiber_t_rec(fiber_t_rec_t0_idx);
                    fiber_itot_rec = fiber_itot_rec(:,fiber_t_rec_t0_idx:end);
                    % fiber_v_rec = fiber_v_rec(:,fiber_t_rec_t0_idx:end);
                end
    
                % ========================= compute SUAP of the fiber =========================
                % LFM           [V/uA]
                % fiber_itot    [mA/cm^2]
                % fiber_area    [um^2]
                % SUAP          [uV]
                SUAP{i} = sum(LFM_fiber_padding .* fiber_itot_rec .* fiber_area') * 10; % uV
                SUAP_t{i} = fiber_t_rec; 
                SUAP_dt(i) = fiber_dt;
                [~, SUAP_idx_max(i)] = max(SUAP{i});
                [~, SUAP_idx_min(i)] = min(SUAP{i});
                SUAP_template_type(i) = fiber_template_type;
                SUAP_CV(i) = fiber_CV;
                SUAP_align_init_spike(i) = align_init_spike;
                SUAP_t0(i) = fiber_t_rec_t0;
                SUAP_t0_idx(i) = fiber_t_rec_t0_idx;

                % ============== DEBUGGING PLOTS - travelling spike vs SUAP ===============
                % figure;
                % pos_z_as = 0; % CHANGE THIS!
                % [~, idx] = min(abs(pos_z_as-locs_fiber_pos_z(locs_fiber_idcs_ranvier)));
                % subplot(2,1,1);
                % plot(fiber_t_rec, fiber_v_rec(locs_fiber_idcs_ranvier(1),:));
                % hold on;
                % plot(fiber_t_rec, fiber_v_rec(locs_fiber_idcs_ranvier(idx),:));
                % yyaxis right;
                % plot(fiber_t_rec, SUAP{i});
                % legend('Ranvier at -nerve length/2', 'Ranvier closest to AS', 'SUAP');
                % subplot(2,1,2);
                % plot(fiber_t_rec, fiber_itot_rec(locs_fiber_idcs_ranvier(1),:));
                % hold on;
                % plot(fiber_t_rec, fiber_itot_rec(locs_fiber_idcs_ranvier(idx),:));
                % linkaxes(findall(gcf,'type','axes'),'x');
                % ==========================================================================
    
            end
            
            % store results
            SUAP_data.SUAP = SUAP;
            SUAP_data.SUAP_t = SUAP_t;
            SUAP_data.SUAP_dt = SUAP_dt;
            SUAP_data.SUAP_idx_max = SUAP_idx_max;
            SUAP_data.SUAP_idx_min = SUAP_idx_min;
            SUAP_data.SUAP_template_type = SUAP_template_type;
            SUAP_data.SUAP_CV = SUAP_CV;
            SUAP_data.SUAP_align_init_spike = SUAP_align_init_spike;
            SUAP_data.SUAP_t0 = SUAP_t0;
            SUAP_data.SUAP_t0_idx = SUAP_t0_idx;
        end

        function [topo_func_sel, topo_func_sel_idcs] = get_slice_topo_func(obj, id_fibers)
            % Extract a given set of fibers (defined by their ids) from topo_func 
            % 
            % [topo_func_sel, topo_func_sel_idcs] = get_slice_topo_func(id_fibers)

            topo_func_sel_idcs = zeros(size(id_fibers));
            for i=1:length(id_fibers)
                topo_func_sel_idcs(i) = find(id_fibers(i) == obj.topo_func(:,1));
            end
            topo_func_sel = obj.topo_func(topo_func_sel_idcs,:);
        end

    end
end