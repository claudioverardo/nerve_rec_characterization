%% SIMULATE DATA USED IN FIGURE 6C
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber_cuff_coverage");
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
% opts_i{11}    y_el        [um]    y-position of the electrode (only for TIME/uNG)
% opts_i{12}    H_AS        [mm]    height of the active site (size along the cuff circumference)
% opts_i{13}    W_AS        [mm]    width of the active site (size along the longitudinal direction)
% ---------------------------------------------------------------------------------------------------------------


% uncomment the code blocks below to create/simulate the desired subset of hybrid models
% =============================================================================================================================
% EXTRANEURAL (cuff)
R_nerve = 1e-3; % m
H_AS = [[0.02 0.1 0.2 0.5 1.0]*2*pi*R_nerve*1e+3];
W_AS = [ 1.0  1.0 1.0 1.0 1.0                   ];
d_as_fascs = [20];
d_fasc_fib = 200;
% myelinated fibers
d_fibers = [2 4 6 8 10 12 14];
dz = 5; % um
z_shift = 0;
for i=1:length(d_as_fascs)
    d_as_fasc = d_as_fascs(i);
    for j=1:length(H_AS)
        opts{end+1} = {"Cortec_Cuff"    2      4       200     -800+d_as_fasc       nan         "ascent_mrg"         d_fibers         -1000+d_fasc_fib+d_as_fasc       z_shift     nan   H_AS(j)   W_AS(j)};
    end
end
% unmyelinated fibers
d_fibers_unmyel = [0.5 1.0 1.5];
z_shift = 0;
for i=1:length(d_as_fascs)
    d_as_fasc = d_as_fascs(i);
    for j=1:length(H_AS)
        opts{end+1} = {"Cortec_Cuff"    2      4       200     -800+d_as_fasc       nan         "ascent_sundt15"     d_fibers_unmyel  -1000+d_fasc_fib+d_as_fasc       z_shift     nan   H_AS(j)   W_AS(j)};
    end
end
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

        hm_params.electrodes.radius = R_nerve;
        hm_params.electrodes.radius_nerve = R_nerve;

        hm_params.electrodes.length = 15e-3;
        hm_params.electrodes.thick = 1e-3;
        hm_params.electrodes.n_as = 3;
        hm_params.electrodes.length_as = opts_i{12}*1e-3;
        hm_params.electrodes.depth_as = 20e-6;
        hm_params.electrodes.width_as = opts_i{13}*1e-3;
        hm_params.electrodes.pitch_as = 4e-3;

        hm_params.electrodes.z = 0;
        hm_params.electrodes.theta_z = 90;

        % ================================================================================= 

        % mphstart
        hm_model = HMModel(hm_params);
        hm_model = hm_model.FEMModel_build();

        hm_model.export_node_locs('log_step', 100, 'node_length_unmyel', 50/6/10);
        hm_model.export_LFM([1 1 1], [1 2 3], 'mesh_size', opts_i{2});
        hm_model.export_SUAP_template([1 1 1], [1 2 3], dir_templates, 'log_step', 1);

    end
end