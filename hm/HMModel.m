classdef HMModel
% HMModel Fundamental class for the development of hybrid models (HMs).
% 
% HMModel Properties:
%     nerve                             - definition of the nerve
%     implant                           - definition of the implant(s)
%     saline                            - definition of the saline external medium
%     fem_model                         - handle of the interface with Comsol (FEM)
%     dir_model                         - directory where the HM is stored
%   
% HMModel Methods:
% 
%     ------------ HM model initialization -------------
%     HMModel                           - constructor
% 
%     ----------- FEM (interface w/ Comsol) ------------
%     FEMModel_build                    - build the Comsol model and store its mph file
%     FEMModel_init                     - initialize the Comsol model
%     FEMModel_generate_materials       - generate materials in the Comsol model
%     FEMModel_generate_geometry        - generate geometry in the Comsol model
%     FEMModel_assign_materials         - assign materials to geometry domains in the Comsol model
%     FEMModel_remove_elec_from_nerve   - simplify geometry definition (nerve-electrode binary operations) in the Comsol model
%     FEMModel_generate_mesh            - generate mesh in the Comsol model
%     FEMModel_solve_study              - solve the volume conductor problem in the Comsol model
%     FEMModel_get_implants_selections  - get the tags of the selections of the implants subdomains in the Comsol model
%     FEMModel_load                     - load the Comsol model from the stored .mph file 
%     FEMModel_clear                    - remove the Comsol model and its Matlab handle
% 
%     ----------- HM model characterization ------------
%     export_node_locs                  - compute and export the node (compartment) locations of the fibers
%     export_LFM                        - compute and export the lead-field matrices (LFMs) of the active sites
%     export_SUAP_template              - compute and export the single-unit action potentials (SUAPs) recorded by the active sites
% 
%     ----------------- HM model utils -----------------
%     load_locs                         - load the node (compartment) locations of the fibers
%     load_locs_meta                    - load the metadata related to the node (compartment) locations of the fibers
%     load_LFM                          - load the LFMs of the active sites
%     load_SUAP_template                - load the SUAPs of the active sites
%     get_all_elec_site_ids             - get the identifiers of all the electrodes and active sites
%     elec2implant_id                   - convert from unique electrode identifier to implant+electrode identifier 

    properties
        % nerve - definition of the nerve (1x1 PeripheralNerve)
        nerve

        % implant - definition of the implant(s) (1xN Implant)
        implant

        % saline - definition of the saline external medium (1x1 Saline)
        saline
        
        % fem_model - handle of the interface with Comsol (FEM) (1x1 ModelClient -> Comsol LiveLink class)
        fem_model

        % dir_model - directory where the HM is stored (1x1 string)
        dir_model
    end
    
    methods
    
        function obj = HMModel(params)
            % Constructor of the class.
            % Instantiate properties of the HMModel object from a configuration struct.
            % 
            % obj = HMModel(params)
            % 
            % PARAMETERS
            % ----------
            % params (1x1 struct)
            %     struct defining the parameters of the hybrid model
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object

            obj.dir_model = params.dir_model;
            if ~isfolder(obj.dir_model)
                mkdir(obj.dir_model);
            end
            
            % ======================== SALINE ======================= 
            obj.saline = Saline(params.saline);

            % ======================== NERVE ======================== 
            obj.nerve = PeripheralNerve(params);
            
            % ================ IMPLANTS + ELECTRODES ================ 
            n_implants = length(params.implant_type);

            if n_implants == 0
                params.electrodes = [];
            end

            % legacy
            params.n_elecs = ones(size(params.implant_type));

            if length(params.electrodes) ~= sum(params.n_elecs)
                error("electrode parameters must be defined for each electrode");
            end
            
            % create implants
            elecs_count = 0;
            obj.implant = Implant.empty(n_implants,0);
            for i=1:n_implants
                id_elecs = elecs_count + (1:params.n_elecs(i));
                params_implant = struct('implant_type', params.implant_type(i), 'n_elecs', params.n_elecs(i), 'id_elecs', id_elecs, 'electrodes', params.electrodes(id_elecs));
                obj.implant(i) = Implant(params_implant);
                elecs_count = elecs_count + params.n_elecs(i);
            end
            
            % ============= STORE hm_params + hm_model ============== 
            hm_params = params;
            save(fullfile(obj.dir_model,"hm_params.mat"), 'hm_params');
            hm_model = obj;
            save(fullfile(obj.dir_model,"hm_model.mat"), 'hm_model');
        end

        function obj = FEMModel_build(obj)
            % Build the Comsol model from the Matlab objects and store it in a .mph file.
            % It calls the following methods:
            %     1) FEMModel_init()
            %     2) FEMModel_generate_materials()
            %     3) FEMModel_generate_geometry()
            %     4) FEMModel_assign_materials()
            %     5) FEMModel_remove_elec_from_nerve()
            % 
            % obj = FEMModel_build()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            
            disp("======== run fem_model build ========");
            disp(obj.dir_model);

            obj = obj.FEMModel_init();
            obj = obj.FEMModel_generate_materials();
            obj = obj.FEMModel_generate_geometry();
            obj = obj.FEMModel_assign_materials();
            obj = obj.FEMModel_remove_elec_from_nerve();

            mphsave(obj.fem_model.tag, fullfile(obj.dir_model,"fem_model.mph"));
        end

        function obj = FEMModel_init(obj)
            % Instantiate a blank Comsol model and add:
            %     - a 3D component
            %     - the electric currents physics
            %     - a stationary study
            %     - an automatically generated mesh
            % 
            % obj = FEMModel_init()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)

            import com.comsol.model.*
            import com.comsol.model.util.*
            ModelUtil.showProgress(true);
            obj.fem_model = ModelUtil.create('Model');
            obj.fem_model.geom.create('geom1', 3);
            obj.fem_model.physics.create('ec', 'ConductiveMedia', 'geom1');
            obj.fem_model.study.create('std1');
            obj.fem_model.study('std1').create('stat', 'Stationary');
            obj.fem_model.study('std1').feature('stat').activate('ec', true);

            obj.fem_model.sol.create('sol1');
            obj.fem_model.sol('sol1').study('std1');
            obj.fem_model.sol('sol1').create('st1', 'StudyStep');
            obj.fem_model.sol('sol1').create('v1', 'Variables');
            obj.fem_model.sol('sol1').create('s1', 'Stationary');
            obj.fem_model.sol('sol1').feature('s1').create('i1', 'Iterative');
            obj.fem_model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
            obj.fem_model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
            obj.fem_model.sol('sol1').feature('s1').set('stol', '1e-6');

            obj.fem_model.mesh.create('mesh', 'geom1');
            obj.fem_model.mesh('mesh').autoMeshSize(3); % default, the mesh size can be modified with the routines to run simulations
        end
        
        function obj = FEMModel_generate_materials(obj)
            % Add materials in the Comsol model.
            % 
            % obj = FEMModel_generate_materials()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            
            endo = obj.fem_model.material.create('endo');
            endo.name('Endoneurium');
            endo.propertyGroup('def').set('relpermittivity', 80); % water
            endo.propertyGroup('def').set('electricconductivity',{'0.083', '0', '0', '0', '0.083', '0', '0', '0', '0.571'});
            
            peri = obj.fem_model.material.create('peri');
            peri.name('Perineurium');
            peri.propertyGroup('def').set('relpermittivity', 80);
            peri.propertyGroup('def').set('electricconductivity',0.0009);
            
            epi = obj.fem_model.material.create('epi');
            epi.name('Epineurium');
            epi.propertyGroup('def').set('relpermittivity', 80);
            epi.propertyGroup('def').set('electricconductivity',0.083);
            
            sal = obj.fem_model.material.create('sal');
            sal.name('Saline');
            sal.propertyGroup('def').set('relpermittivity', 80);
            sal.propertyGroup('def').set('electricconductivity',2);
            
            plat = obj.fem_model.material.create('as');
            plat.name('Sites');
            plat.propertyGroup('def').set('relpermittivity', 1);
            plat.propertyGroup('def').set('electricconductivity', 1e10);
            
            poly = obj.fem_model.material.create('elec');
            poly.name('Elecshaft');
            poly.propertyGroup('def').set('relpermittivity', 1);
            poly.propertyGroup('def').set('electricconductivity',1e-10);
            
            tungsten = obj.fem_model.material.create('tungsten');
            tungsten.name('Tungsten');
            tungsten.propertyGroup('def').set('relpermittivity', 1);
            tungsten.propertyGroup('def').set('electricconductivity',1.79e+7);
            
            encapsulation = obj.fem_model.material.create('encapsulation');
            encapsulation.name('Encapsulation');
            encapsulation.propertyGroup('def').set('relpermittivity', 80);
            encapsulation.propertyGroup('def').set('electricconductivity', 0.159);

        end
        
        function obj = FEMModel_generate_geometry(obj)
            % Add saline, nerve, and implant geometries in the Comsol model.
            % 
            % obj = FEMModel_generate_geometry()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            
            disp('## Starts adding outernerve geometry to the model ##')
            obj.fem_model = obj.saline.add_to_model(obj.fem_model, obj.nerve);
            disp('## Ends adding outernerve geometry to the model ##')
            disp('## Starts adding nerve geometry to the model ##')
            obj.fem_model = obj.nerve.add_to_model(obj.fem_model);
            disp('## Ends adding nerve geometry to the model ##')
            disp('## Starts adding implant(s) geometry to the model ##')
            for i = 1:length(obj.implant)
                obj.fem_model = obj.implant(i).add_to_model(obj.fem_model);
            end
            disp('## Ends adding implant geometry to the model ##')
        end
        
        function obj = FEMModel_assign_materials(obj)
            % Assign materials to the geometrical domains in the Comsol model.
            % 
            % obj = FEMModel_assign_materials()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            
            % ======================== SALINE ======================= 
            obj.saline.assign_materials(obj.fem_model);

            % ======================== NERVE ======================== 
            obj.nerve.assign_materials(obj.fem_model);

            % ================ IMPLANTS + ELECTRODES ================ 
            % assign materials shafts / active sites of electrodes
            for i = 1:length(obj.implant)
                obj.implant(i).assign_materials(obj.fem_model);
            end
            
        end
        
        function obj = FEMModel_remove_elec_from_nerve(obj)
            % Perfom binary operations to simplify the geometry representation in the Comsol model.
            % Specifically, subtract electrode domains from nerve + saline domains.
            %
            % obj = FEMModel_remove_elec_from_nerve()
            % 
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            
            disp('## Starts simplifying geometry ##')
                
            geom1 = obj.fem_model.geom('geom1');
            
            [elec_sel, as_sel, encaps_sel] = FEMModel_get_implants_selections(obj);
            if ~isempty(elec_sel)
                elec_sel = [elec_sel{:}];
                elec_sel = [elec_sel{:}];
            end
            if ~isempty(as_sel)
                as_sel = [as_sel{:}];
                as_sel = [as_sel{:}];
            end
            if ~isempty(encaps_sel)
                encaps_sel = [encaps_sel{:}];
                encaps_sel = [encaps_sel{:}];
            end
            elec_as_sel = [elec_sel, as_sel, encaps_sel];
            
            saline = geom1.create('saline', 'Difference');
            saline.selection('input').set('salgeom');
            saline.selection('input2').set(elec_as_sel);
            saline.set("keepsubtract", true);
            saline.set('createselection', 'on');

            epineurium = geom1.create('epineurium', 'Difference');
            epineurium.selection('input').set('epigeom');
            epineurium.selection('input2').set(elec_as_sel);
            epineurium.set("keepsubtract", true);
            epineurium.set('createselection', 'on');

            endoneurium = geom1.create('endoneurium', 'Difference');
            endoneurium.selection('input').set('endogeom');
            endoneurium.selection('input2').set(elec_as_sel);
            endoneurium.set("keepsubtract", true);
            endoneurium.set('createselection', 'on');
            
            geom1.run();

            disp('## Ends simplifying geometry ##')
        end

        function obj = FEMModel_generate_mesh(obj, varargin)
            % Generate the mesh in the Comsol model.
            % 
            % obj = FEMModel_generate_mesh(varargin)
            % 
            % NAME, VALUE ARGS
            % ----------------
            % mesh_type (string)
            %     'auto' (default), 'custom' (TODO)
            % mesh_size (int)
            %     set the element size of the physics-controlled mesh 
            %     1: Extremely fine
            %     2: Extra fine
            %     3: Finer (default)
            %     4: Fine
            %     5: Normal
            %     6: Coarse
            %     7: Coarser
            %     8: Extra coarse
            %     9: Extremely coarse
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)

            % input parser
            p = inputParser;
            addParameter(p, 'mesh_type', 'auto');
            addParameter(p, 'mesh_size', 3);
            p.KeepUnmatched = true;
            parse(p, varargin{:});
            
            % parse function parameters
            mesh_type = p.Results.mesh_type;
            mesh_size = p.Results.mesh_size;

            if ~strcmp(mesh_type, 'auto')
                error("mesh_type not supported")
            end

            obj.fem_model.mesh('mesh').autoMeshSize(mesh_size);
            obj.fem_model.mesh('mesh').run;
        end
        
        function [obj, id] = FEMModel_solve_study(obj, implant_id, implant_elec_ids, site_ids, currs, fem_basefilenames, coord_filenames, id, varargin)
            % Solve the volume conductor problem (stationary study) in the Comsol model.
            % 
            % [obj, id] = FEMModel_solve_study(implant_id, implant_elec_ids, site_ids, currs, fem_basefilenames, coord_filenames, id)
            %
            % PARAMETERS
            % ----------
            % implant_id (n_actsites x 1 int array)
            %     identifier of the implant for each active site
            % implant_elec_ids (n_actsites x 1 int array)
            %     identifier of the electrode inside the implant for each active site
            % site_ids (n_actsites x 1 int array)
            %     identifier for each active site
            % currs (n_actsites x 1 double array)
            %     stimulating current for each active site
            % fem_basefilenames
            %     basefilenames of the files where to export the sampled FEM solutions
            % coord_filenames
            %     basefilenames of the files defining the sampling points of the FEM solutions
            % id (int)
            %     identifier for the current simulation
            % 
            % NAME, VALUE ARGS
            % ----------------
            % see FEMModel_generate_mesh()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            % id (int)
            %     identifier for the next stimulation
            
            % generate mesh 
            obj.FEMModel_generate_mesh(varargin{:});
            % obj.fem_model.mesh('mesh').autoMeshSize(mesh_size);
            % obj.fem_model.mesh('mesh').run;

            % set all simulating currents
            for i = 1:length(implant_elec_ids)
                obj.fem_model = obj.implant(implant_id).electrodes(implant_elec_ids(i)).set_active_site_current(obj.fem_model, site_ids(i), currs(i));
            end
            
            % run stationary study with the default solver
            % obj.fem_model.study('std1').run;

            % run stationary study with the custom solver (see FEMModel_init())
            obj.fem_model.sol('sol1').runAll;

            % export electric potential in mesh points or interpolated to
            % specified locations
            if isempty(coord_filenames)
                dataref = mpheval(obj.fem_model,{'V'});
                save(fem_basefilenames, 'dataref');
            else
                if ~iscell(coord_filenames) && ~isstring(coord_filenames)
                    coord_filenames = {coord_filenames};
                    fem_basefilenames = {fem_basefilenames};
                end
                for i = 1:length(coord_filenames)
                    data = obj.fem_model.result.export.create(strcat('data_', num2str(i), '_', int2str(id)), 'Data');
                    data.setIndex('expr', 'V', 0);
                    data.set('filename', strcat(fem_basefilenames{i}, '.txt'));
                    data.set('location', 'file');
                    data.set('coordfilename', strcat(coord_filenames{i}, '.txt'));
                    data.set('fullprec', true);
                    data.run;
                    obj.fem_model.result.export.remove(data.tag);
                end
            end
            id = id + 1;
            disp('Potential distribution computed and exported')

            % reset all stimulating currents to zero (for next simulation)
            for i = 1:length(implant_elec_ids)
                obj.fem_model = obj.implant(implant_id).electrodes(implant_elec_ids(i)).set_active_site_current(obj.fem_model, site_ids(i), 0);
            end
        end

        function [elec_sel, as_sel, encaps_sel] = FEMModel_get_implants_selections(obj)
            % Retrieve the tags of the selections to access the subdomains of the implants. 
            % 
            % [elec_sel, as_sel, encaps_sel] = FEMModel_get_implants_selections()
            %
            % RETURNS
            % -------
            % elec_sel (cell array)
            %     tags of the selections for the electrode shafts
            %     elec_sel{i}{j} refers to the i-th implant and its j-th electrode
            % as_sel (cell array)
            %     tags of the selections for the active sites
            %     as_sel{i}{j} refers to the i-th implant and its j-th electrode
            % encaps_sel (cell array)
            %     tags of the selections for the encapsulation tissue
            %     encaps_sel{i}{j} refers to the i-th implant and its j-th electrode

            elec_sel = {};
            as_sel = {};
            encaps_sel = {};
            for i = 1:length(obj.implant)
                for j = 1:obj.implant(i).n_elecs
                    [elec_sel{i}{j}, as_sel{i}{j}, encaps_sel{i}{j}] = obj.implant(i).electrodes(j).get_selections();
                end
            end
        end

        function obj = FEMModel_load(obj)
            % Load the Comsol model from the saved .mph file.
            % 
            % obj = FEMModel_load()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)

            obj.fem_model = mphload(fullfile(obj.dir_model,'fem_model.mph'));
        end

        function obj = FEMModel_clear(obj)
            % Remove the Comsol model (and set obj.fem_model = []).
            % 
            % obj = FEMModel_clear()
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)

            import com.comsol.model.*
            import com.comsol.model.util.*
            ModelUtil.remove(obj.fem_model.tag);
            obj.fem_model = [];
        end

        function obj = export_node_locs(obj, varargin)
            % Generate 3D locations of fiber nodes (compartments) from the functional topography. 
            % Store them and their metadata in files in the hm_model directory.
            % 
            % NOTE: wrapper to obj.nerve.get_node_locs() and obj.nerve.fiber_topography.get_node_locs()
            %       to automatically handle the length of the nerve and all the fibers in the nerve.
            % 
            % obj = export_node_locs(varargin)
            % 
            % NAME, VALUE ARGS
            % ----------------
            % node_length_unmyel
            %     spatial discretization step for unmyelinated fibers (50/6 um by default)
            % log_step (int)
            %     control amount of log information in the command window
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated COMSOL model)
            %
            % EXPORT FILE 1: locs.txt
            % -----------------------
            % 3D locations of the nodes of the fibers. For each fiber, the node locations are generated
            % along the longitudinal direction of the nerve from -obj.nerve.length/2 to +obj.nerve.length/2.
            % Then, the node locations of all the fibers are vertically concatenated.
            %     locs(:,1)   (n_locs x 1 double)    x-location of fiber nodes [m]
            %     locs(:,2)   (n_locs x 1 double)    y-location of fiber nodes [m]
            %     locs(:,3)   (n_locs x 1 double)    z-location of fiber nodes [m]
            % 
            % NOTE: n_locs = sum_i^n_fibers n_nodes(i), where n_nodes(i) is the number of nodes of 
            %       the i-th fiber. n_nodes(i) is not a parameter directly controlled but it emerges
            %       indirectly from the spatial discretization rule of the fiber, the longitudinal 
            %       shift of the fiber, and the length of the nerve.
            %
            % EXPORT FILE 2: locs_meta.mat
            % ----------------------------
            % metadata related to the generated 3D locations of fiber nodes
            %     locs_id_fiber               (n_locs x 1 double)
            %         id of the fiber associated to each node location
            %     locs_mask_ranvier           (n_locs x 1 bool)
            %         boolean mask denoting which node locations are Ranvier nodes
            %     n_nodes                     (n_fibers x 1 double)
            %         number of node locations generated for each fiber
            %     z_shift_m                   (n_fibers x 1 double)
            %         longitudinal shift of each fiber [m]
            %     node_length_unmyel          (1 x 1 double)
            %         spatial discretization step used for unmyelinated fibers
            %     n_internodes_myel           (n_fibers x 1 double)
            %         number of internodes to be instantiated for each myelinated fiber in Neuron simulations
            %         (if the i-th fiber is unmyelinated, n_internodes_myel(i) = 0).
            %         WARNING: in NEURON simulations, the instantiated fibers can be longer than the nerve.
            %                  Thus, the number of nodes generated for the i-th fiber can be > n_nodes(i).
            %                  Do not use this array in Matlab routines.
            %     n_nodes_unmyel              (n_fibers x 1 double)
            %         number of nodes to be instantiated for each unmyelinated fiber in Neuron simulations
            %         (if the i-th fiber is myelinated, n_nodes_unmyel(i) = 0).
            %         WARNING: in NEURON simulations, the instantiated fibers can be longer than the nerve.
            %                  Thus, the number of nodes generated for the i-th fiber can be > n_nodes(i).
            %                  Do not use this array in Matlab routines.
            %     n_internodes_myel_step      (1 x 1 double)
            %         DEPRECATED - do not use
            %     n_nodes_unmyel_step         (1 x 1 double)
            %         DEPRECATED - do not use
            %     locs_idcs_inside_nerve      (n_locs x 1 double)
            %         DEPRECATED - do not use
            
            disp("========== run locs export ==========");
            disp(fullfile(obj.dir_model,'locs.txt'));

            % ============================= generate locs =============================
            [node_locs, node_locs_metadata] = obj.nerve.get_node_locs([],varargin{:});
    
            % ============================== store locs ===============================
            writematrix(node_locs, fullfile(obj.dir_model,"locs.txt"));
        
            % =================== store metadata of locs generation ===================
            save(fullfile(obj.dir_model,"locs_meta.mat"),'-struct','node_locs_metadata');
        
        end

        function obj = export_LFM(obj, elec_ids, site_ids, varargin)
            % Compute the lead-field matrix (LFM) for a set of active sites. For each site, the
            % volume conductor problem is solved in Comsol and the results are sampled a the 
            % fiber nodes. The solution of the volume conductor model is not influenced by the
            % presence of fibers. Store results in files in the hm_model directory.
            % 
            % NOTE: wrapper to FEMModel_solve_study() to automatically handle the active sites in
            %       the model and the correspondent inport/export of the locs/LFM_#elec_#as files.
            %       Further, it enables the use of a temporary dir to circumvent the limitations
            %       of Comsol/Windows on the path length.
            % 
            % obj = export_LFM(elec_ids, site_ids, varargin)
            % 
            % PARAMETERS
            % ----------
            % elec_ids (n_actsites x 1 int array)
            %     electrode identifier for each active site
            % site_ids (n_actsites x 1 int array)
            %     identifier for each active site
            % 
            % NAME, VALUE ARGS
            % ----------------
            % mesh_size (int)
            %     set the element size of the physics-controlled mesh (1-9, set to 3 by default)
            %     see FEMModel_generate_mesh()
            % dir_temp (string)
            %     temporary directory for LFM export (to circumvent Comsol limit on path lengths)
            %
            % RETURNS
            % -------
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            %
            % EXPORT FILES: LFM_[elec_ids(i)]_[site_ids(i)].txt
            % -------------------------------------------------
            % The i-th file stores the LFM of the i-th active site provided in input, as follows:
            %     locs(:,1)   (n_locs x 1 double)    x-location of fiber nodes [m]
            %     locs(:,2)   (n_locs x 1 double)    y-location of fiber nodes [m]
            %     locs(:,3)   (n_locs x 1 double)    z-location of fiber nodes [m]
            %     locs(:,4)   (n_locs x 1 double)    LFM @ the fiber nodes [V/uA]
            
            % input parser
            p = inputParser;
            addParameter(p, 'mesh_size', 3);
            addParameter(p, 'dir_temp', '');
            p.KeepUnmatched = true;
            parse(p, varargin{:});
            
            % parse function parameters
            mesh_size = p.Results.mesh_size;
            dir_temp = p.Results.dir_temp;
            % args_unmatched = p.Unmatched;

            disp("=========== run LFM export ==========");

            locs_filename = "locs";
            LFM_filename = "LFM";

            [implant_ids, implant_elec_ids] = obj.elec2implant_id(elec_ids);
            n_sites = length(site_ids);

            % ========================== compute LFM (w/ FEM) ==========================
            % NOTE: require obj.fem_model loaded and Comsol running
            if isempty(obj.fem_model)
                error("fem_model not found, use FEMModel_build() or FEMModel_load()");
            end

            if strcmp(dir_temp,'')
                USE_TEMP_DIR = false;
            else
                USE_TEMP_DIR = true;
            end
        
            % workaround to handle very long paths in Comsol data export
            if USE_TEMP_DIR
                [~, name_hm, ~] = fileparts(obj.dir_model);

                if use_locs
                    locs_filename_temp = name_hm+"_"+locs_filename;
                    copyfile(fullfile(obj.dir_model,locs_filename+".txt"), fullfile(dir_temp,locs_filename_temp+".txt"));
                end

                LFM_filename_temp = name_hm+"_"+LFM_filename;
            end

            for n = 1:n_sites
                k = implant_ids(n);         % idx of the implant
                i = implant_elec_ids(n);    % idx of the electrode inside the implant
                j = site_ids(n);            % idx of the site inside the electrode

                disp(strcat('# Studying-', obj.implant(k).implant_type, ' electrode id.', num2str(obj.implant(k).electrodes(i).id), ', site n.', num2str(j)))
                LFM_filename_out = strcat(LFM_filename, '_', num2str(obj.implant(k).electrodes(i).id) , '_', num2str(j));  
                disp(fullfile(obj.dir_model,LFM_filename_out+".txt")); 

                id = 1;
                
                % workaround to handle very long paths in Comsol data export
                if USE_TEMP_DIR
                    [obj, id] = obj.FEMModel_solve_study(k, i, j, 1e-6, fullfile(dir_temp,LFM_filename_temp), fullfile(dir_temp,locs_filename_temp), id, 'mesh_type', 'auto', 'mesh_size', mesh_size);
                    movefile(fullfile(dir_temp,LFM_filename_temp+".txt"), fullfile(obj.dir_model,LFM_filename_out+".txt"));
                else
                    [obj, id] = obj.FEMModel_solve_study(k, i, j, 1e-6, fullfile(obj.dir_model,LFM_filename_out), fullfile(obj.dir_model,locs_filename), id, 'mesh_type', 'auto', 'mesh_size', mesh_size);
                end
            end
            % workaround to handle very long paths in Comsol data export
            if USE_TEMP_DIR && use_locs
                delete(fullfile(dir_temp,locs_filename_temp+".txt"));
            end

        end

        function obj = export_SUAP_template(obj, elec_ids, site_ids, dir_fiber_templates, varargin)
            % Compute the single-unit action potentials (SUAPs) of the fibers for a set 
            % of active sites. For each site, the SUAPs are obtained by multiplying the
            % templates of fiber excitations by the LFM of the site. The results are then 
            % stored files in the hm_model directory.
            % 
            % Use excitation_template_fibers_parallel.py to generate the templates.
            % An excitation template consists of the time evolution of membrane potential  
            % (v) and transmembrane current (itot) during the propagation along the fiber
            % of a spike evoked by intracellular stimulation (IClamp, synapses, etc).
            % 
            % Two types of excitation templates can be generated/handled:
            % - "internode-CV": v/itot are stored only at the nodes of an internodal 
            %                   region, along with the conduction velocity (CV). 
            %                   The v/itot of all the nodes are obtained as time-shifted
            %                   versions according to the CV.
            % - "full-fiber":   v/itot are stored at all the nodes of the fiber. 
            %                   The fiber in the excitation template is assumed longer
            %                   than the nerve.
            % 
            % The code assumes an excitation template to be available for each fiber 
            % in the functional topography (discrete set of diameters admitted).
            % 
            % NOTE: this function wrappers obj.nerve.fiber_topography.get_SUAP_template()
            %       to automatically handle all the fibers and specific active sites in the
            %       model, with the correspondent inport/export of the locs, locs_meta,
            %       LFM_#elec_#as, SUAP_#elec_#as files.
            % 
            % obj = export_SUAP_template(elec_ids, site_ids, dir_fiber_templates, varargin)
            %             
            % PARAMETERS
            % ----------
            % elec_ids (n_actsites x 1 int array)
            %     electrode identifier for each active site
            % site_ids (n_actsites x 1 int array)
            %     identifier for each active site
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
            % obj (1x1 HMModel)
            %     the modified caller object (w/ updated Comsol model)
            %
            % EXPORT FILES: SUAP_[elec_ids(i)]_[site_ids(i)].mat
            % --------------------------------------------------
            % The i-th file stores the SUAPs of the i-th active site provided in input, as follows:
            %     SUAP                    (n_fibers x 1 cell array)     [uV]
            %     SUAP_t                  (n_fibers x 1 cell array)     [ms]
            %     SUAP_dt                 (n_fibers x 1 double)         [ms]
            %     SUAP_idx_max            (n_fibers x 1 double)         [-]
            %     SUAP_idx_min            (n_fibers x 1 double)         [-]
            %     SUAP_template_type      (n_fibers x 1 double)         [-]
            %     SUAP_CV                 (n_fibers x 1 double)         [m/s]
            %     SUAP_align_init_spike   (1 x 1 bool]                  [-]
            %     SUAP_t0                 (n_fibers x 1 double)         [ms]
            %     SUAP_t0_idx             (n_fibers x 1 double)         [ms]
            
            % input parser
            p = inputParser;
            addParameter(p, 'v_thresh', -40);
            addParameter(p, 'log_step', 50);
            parse(p, varargin{:});
            
            % parse function parameters
            v_thresh = p.Results.v_thresh;
            log_step = p.Results.log_step;
            
            disp("===== run template SUAP export ======");

            id_fibers = obj.nerve.fiber_topography.topo_func(:,1);
            load(fullfile(obj.dir_model,"locs_meta.mat"), 'z_shift_m', 'locs_id_fiber', 'locs_mask_ranvier');
            
            n_sites = length(elec_ids);
            for n = 1:n_sites
                
                % ============================= compute SUAP ==============================
                LFM = readmatrix(fullfile(obj.dir_model,"LFM_"+elec_ids(n)+"_"+site_ids(n)+".txt"));
                SUAP_filename = "SUAP_"+elec_ids(n)+"_"+site_ids(n)+".mat";
                SUAP_filepath = fullfile(obj.dir_model,SUAP_filename);
                disp(SUAP_filepath);
                SUAP_data = obj.nerve.fiber_topography.get_SUAP_template( ...
                    id_fibers, obj.nerve.length, LFM(:,1:3), z_shift_m, locs_id_fiber, locs_mask_ranvier, LFM(:,4), dir_fiber_templates, ...
                    'v_thresh', v_thresh, 'log_step', log_step);
                    
                % ============================ store SUAP data ============================
                save(SUAP_filepath,'-struct','SUAP_data','-v7.3');
            end
        
        end

        function locs = load_locs(obj)
            % Return the 3D locations of the fiber nodes
            % 
            % locs = load_locs()

            filename = fullfile(obj.dir_model,"locs.txt");
            if ~isfile(filename)
                error("locs not exported/found: "+filename);
            end
            locs = readmatrix(filename); 
        end

        function locs_meta = load_locs_meta(obj)
            % Return the metadata of the 3D locations of the fiber nodes
            % 
            % locs_meta = load_locs_meta()

            filename = fullfile(obj.dir_model,"locs_meta.mat");
            if ~isfile(filename)
                error("locs_meta not exported/found: "+filename);
            end
            locs_meta = load(filename); 
        end

        function LFM = load_LFM(obj, elec_id, site_id)
            % Return the LFM of an active site
            % 
            % LFM = load_LFM(elec_id, site_id)

            filename = fullfile(obj.dir_model,"LFM_"+elec_id+"_"+site_id+".txt");
            n_head_lines = 9;

            if ~isfile(filename)
                error("LFM not exported/found: "+filename);
            end
            LFM = readmatrix(filename, "NumHeaderLines", n_head_lines); 
        end

        function SUAP = load_SUAP_template(obj, elec_id, site_id)
            % Return the SUAP templates of an active site
            % 
            % SUAP = load_SUAP_template(elec_id, site_id)

            filename = fullfile(obj.dir_model,"SUAP_"+elec_id+"_"+site_id+".mat");

            if ~isfile(filename)
                error("SUAP not exported/found: "+filename);
            end
            SUAP = load(filename); 
        end

        function [elec_ids, site_ids] = get_all_elec_site_ids(obj)
            % Return all the active sites in the model
            % 
            % [elec_ids, site_ids] = get_all_elec_site_ids()
            % 
            % elec_ids => identifiers of the active sites' electrodes
            % site_ids => identifiers of the active sites inside their electrodes
            
            elec_ids = [];
            site_ids = [];

            for k = 1:length(obj.implant)
                for i = 1:obj.implant(k).n_elecs
                    site_ids_aux = 1:obj.implant(k).electrodes(i).n_sites;
                    elec_ids = [elec_ids obj.implant(k).electrodes(i).id*ones(size(site_ids_aux))];
                    site_ids = [site_ids site_ids_aux];
                end
            end
        end

        function [implant_ids, implant_elec_ids] = elec2implant_id(obj, elec_ids)
            % Convert from unique electrode identifier to implant+electrode identifier  
            % 
            % [implant_ids, implant_elec_ids] = elec2implant_id(elec_ids)

            implant_ids = zeros(size(elec_ids));
            implant_elec_ids = zeros(size(elec_ids));

            for n = 1:length(elec_ids)
                for k = 1:length(obj.implant)
                    for i = 1:obj.implant(k).n_elecs
                        if obj.implant(k).electrodes(i).id == elec_ids(n)
                            implant_ids(n) = k;
                            implant_elec_ids(n) = i;
                        end
                    end
                end
            end
        end
    
    end
end