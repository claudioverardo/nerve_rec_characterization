%% SIMULATE DATA USED IN FIGURES 1G, S1
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis");
dir_templates = fullfile(DIR_DATA, "template_fibers");

opts = {};
% opts_i = opts{i,:} defines the i-th hybrid model to be instantiated and simulated
% the model is stored in a subdirectory of dir_experiments named with a MD5 hash string derived from opts_i
% ---------------------------------------------------------------------------------------------------------------
%               parameter   units   description
% ---------------------------------------------------------------------------------------------------------------
% opts_i{1}     electrode           electrode type
% opts_i{2}     mesh                mesh size 
% opts_i{3}     ell         [cm]    extrusion length of the nerve
% opts_i{4}     r_fasc      [um]    radius of the fascicle
% opts_i{5}     y_fasc      [um]    y-position of the fascicle (center)
% opts_i{6}     d_as_fasc   [um]    distance between active site and fascicle (only for extrafascicular TIME/uNG)
% opts_i{7}     mod_fibers          biophysical model of fibers
% opts_i{8}     d_fibers    [um]    diameters of fibers
% opts_i{9}     y_fibers    [um]    y-position of fibers
% opts_i{10}    z-shift             longitudinal shift of fibers
% ---------------------------------------------------------------------------------------------------------------


% uncomment the code blocks below to create/simulate the desired subset of hybrid models
% =============================================================================================================================
% intrafascicular TIME - myelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     20          0.5};
% =============================================================================================================================
% intrafascicular TIME - unmyelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       nan         "ascent_sundt15"     [0.1:0.1:1.5]  20          0.5};
% =============================================================================================================================
% intrafascicular uNG - myelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     20          0.5};
% =============================================================================================================================
% intrafascicular uNG - unmyelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       nan         "ascent_sundt15"     [0.1:0.1:1.5]  20          0.5};
% =============================================================================================================================
% extrafascicular TIME - myelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       10          "ascent_mrg"         [1:0.5:13]     -200+20-10  0.5};
% =============================================================================================================================
% extrafascicular TIME - unmyelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       10          "ascent_sundt15"     [0.1:0.1:1.5]  -200+20-10  0.5};
% =============================================================================================================================
% extrafascicular uNG - myelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       10          "ascent_mrg"         [1:0.5:13]     -200+20-10  0.5};
% =============================================================================================================================
% extrafascicular uNG - unmyelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       10          "ascent_sundt15"     [0.1:0.1:1.5]  -200+20-10  0.5};
% =============================================================================================================================
% extraneural cuff - myelinated fibers
opts{end+1} = {"Cortec_Cuff"   2      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     0           0.5};
% =============================================================================================================================
% extraneural cuff - unmyelinated fibers
opts{end+1} = {"Cortec_Cuff"   2      4       200     0       nan         "ascent_sundt15"     [0.1:0.1:1.5]  0           0.5};
% =============================================================================================================================


% create hybrid models + simulate SUAPs
for i=1:length(opts)
    opts_i = opts{i};
    disp(opts_i);

    hm_params = [];
    hm_params.opts = opts_i;
    hm_params.dir_model = fullfile(dir_experiments, DataHash(opts_i, 'MD5')); 

    if ~(isfile(fullfile(hm_params.dir_model, "hm_model.mat")) && isfile(fullfile(hm_params.dir_model, "SUAP_1_1.mat")))

        % ============================== Bath solution ============================== 
        hm_params.saline.radius_ext = 5 * 1e-3;
        hm_params.saline.topbot = 0;

        % ========================== Nerve structural topography ========================== 
        hm_params.nerve_extrusion_length = opts_i{3}*1e-2;

        R_nerve = 1e-3;
        hm_params.radius_nerve = R_nerve;

        R_fascicle = opts_i{4}*1e-6;
        y_fascicle = opts_i{5}*1e-6;
        hm_params.fascicles = [0 y_fascicle R_fascicle];
        hm_params.perineurium_type = "pelot20_human";

        % ========================== Nerve functional topography ========================== 
        fiber_model_id = get_fiber_model_id(opts_i{7});
        x_fibers = 0;
        y_fibers = opts_i{9}(:)*1e-6;
        diameters = opts_i{8}(:);
        z_shifts = opts_i{10}(:);
        n_fibers = length(y_fibers)*length(diameters)*length(z_shifts);
        id_fibers = (1:n_fibers)';

        % topofunc: id x y mod_id diam z_shift
        topo_func = [zeros(length(y_fibers),1) y_fibers repmat(fiber_model_id,length(y_fibers),1)];
        topo_func = [repmat(topo_func,length(diameters),1) repelem(diameters,size(topo_func,1),1)];
        topo_func = [repmat(topo_func,length(z_shifts),1) repelem(z_shifts,size(topo_func,1),1)];
        topo_func = [id_fibers topo_func];

        hm_params.fibers.topo_func = topo_func;

        % ================================ Electrode ================================ 
        hm_params.implant_type = opts_i{1};

        switch hm_params.implant_type

            case 'Cortec_Cuff'

                hm_params.electrodes.radius = R_nerve;
                hm_params.electrodes.radius_nerve = R_nerve;
                
                hm_params.electrodes.length = 15e-3;
                hm_params.electrodes.thick = 1e-3;
                hm_params.electrodes.n_as = 3;
                hm_params.electrodes.length_as = 2*pi*R_nerve;
                hm_params.electrodes.depth_as = 20e-6;
                hm_params.electrodes.width_as = 1e-3;
                hm_params.electrodes.pitch_as = 4e-3;

                hm_params.electrodes.z = 0;
                hm_params.electrodes.theta_z = 0;

            case 'TIME'

                hm_params.electrodes.l_shaft = 4e-3;
                hm_params.electrodes.w_shaft = 750e-6;
                hm_params.electrodes.h_shaft = 2*9.5e-6;

                hm_params.electrodes.l_cc = 235e-6;
                hm_params.electrodes.n_as = 1;
                hm_params.electrodes.d_as = 80e-6;
                hm_params.electrodes.h_as = 4e-6;

                hm_params.electrodes.x = 0;
                if isnan(opts_i{6})
                    hm_params.electrodes.y = -hm_params.electrodes.h_shaft/2;
                else
                    hm_params.electrodes.y = -hm_params.electrodes.h_shaft/2 - R_fascicle - opts_i{6}*1e-6;
                end
                hm_params.electrodes.z = 0;

                hm_params.electrodes.theta_x = 180; % 1st AS in the y>0 direction
                hm_params.electrodes.theta_y = 0;
                hm_params.electrodes.theta_z = 0;

            case 'uNG'

                hm_params.electrodes.h_shaft = 1.0*hm_params.saline.radius_ext;
                hm_params.electrodes.d_shaft = 200e-6;

                hm_params.electrodes.slope_tip = 5.758;
                hm_params.electrodes.d_rounding_tip = 5e-6;
                hm_params.electrodes.h_uncoated_tip = 30e-6;

                hm_params.electrodes.x = 0;
                if isnan(opts_i{6})
                    hm_params.electrodes.y = 0;
                else
                    hm_params.electrodes.y = -R_fascicle - opts_i{6}*1e-6;
                end
                hm_params.electrodes.z = 0;

                hm_params.electrodes.theta_x = 0;
                hm_params.electrodes.theta_y = 0;
                hm_params.electrodes.theta_z = 0;

        end

        % ================================================================================= 

        % mphstart
        hm_model = HMModel(hm_params);
        hm_model = hm_model.FEMModel_build();

        if strcmp(hm_params.implant_type, 'Cortec_Cuff')
            hm_model.export_node_locs('log_step', 100, 'node_length_unmyel', 50/6/10);
            hm_model.export_LFM([1 1 1], [1 2 3], 'mesh_size', opts_i{2});
            hm_model.export_SUAP_template([1 1 1], [1 2 3], dir_templates, 'log_step', 1);
        else
            hm_model.export_node_locs('log_step', 100);
            hm_model.export_LFM([1], [1], 'mesh_size', opts_i{2});
            hm_model.export_SUAP_template([1], [1], dir_templates, 'log_step', 1);
        end

    end
end
