%% SIMULATE DATA USED IN FIGURES 2-4, 5A-D, 6A, S2-S4
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber");
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
% opts_i{8}     d_fibers    [um]    diameter of fibers
% opts_i{9}     y_fibers    [um]    y-position of fibers
% opts_i{10}    z-shift             longitudinal shifts of fibers
% opts_i{11}    y_el        [um]    y-position of the electrode (only for TIME/uNG)
% ---------------------------------------------------------------------------------------------------------------


% uncomment the code blocks below to create/simulate the desired subset of hybrid models
% =============================================================================================================================
% SPATIAL MAPS - INTRAFASCICULAR (uNG, TIME)
y_els = -[0:2:100]; 
% myelinated fibers
d_fibers = [4 8];
dz = 5; % um
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
    for j=1:length(y_els)
        y_el = y_els(j);
        opts{end+1} = {"uNG"            2      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     y_el};
        opts{end+1} = {"TIME"           1      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     y_el};
    end
end
% unmyelinated fibers
d_fibers = [1];
z_shift = [0; 0.5];
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    for j=1:length(y_els)
        y_el = y_els(j);
        opts{end+1} = {"uNG"            2      4       200     0       nan         "ascent_sundt15"         d_fiber         0          z_shift     y_el};
        opts{end+1} = {"TIME"           1      4       200     0       nan         "ascent_sundt15"         d_fiber         0          z_shift     y_el};
    end
end
% =============================================================================================================================
% SPATIAL MAPS - EXTRAFASCICULAR (uNG, TIME)
y_els = -190-[0:2:110];
% myelinated fibers
d_fibers = [4 8];
dz = 5; % um
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
    for j=1:length(y_els)
        y_el = y_els(j);
        % opts_i{6} = 10 not used !!!
        opts{end+1} = {"uNG"         2      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     y_el};
        opts{end+1} = {"TIME"        1      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     y_el};
    end
end
% unmyelinated fibers
d_fibers = [1];
z_shift = [0; 0.5];
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    for j=1:length(y_els)
        y_el = y_els(j);
        % opts_i{6} = 10 not used !!!
        opts{end+1} = {"uNG"         2      4       200     0       10          "ascent_sundt15"         d_fiber         -190       z_shift     y_el};
        opts{end+1} = {"TIME"        1      4       200     0       10          "ascent_sundt15"         d_fiber         -190       z_shift     y_el};
    end
end
% =============================================================================================================================
% SPATIAL MAPS - EXTRANEURAL (cuff)
d_as_fascs = [0:20:200 300:100:800];
d_fasc_fib = 200;
% myelinated fibers
d_fibers = [4 8];
dz = 5; % um
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
    for j=1:length(d_as_fascs)
        d_as_fasc = d_as_fascs(j);
        opts{end+1} = {"Cortec_Cuff"        2      4       200     -800+d_as_fasc       nan         "ascent_mrg"         d_fiber         -1000+d_fasc_fib+d_as_fasc       z_shift     nan};
    end
end
% unmyelinated fibers
d_fibers = [1];
z_shift = [0; 0.5];
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    for j=1:length(d_as_fascs)
        d_as_fasc = d_as_fascs(j);
        opts{end+1} = {"Cortec_Cuff"        2      4       200     -800+d_as_fasc       nan         "ascent_sundt15"    d_fiber         -1000+d_fasc_fib+d_as_fasc       z_shift     nan};
    end
end
% =============================================================================================================================
% OVERLAP WINDOWS - INTRAFASCICULAR (uNG, TIME)
y_els = -[2 50];
% myelinated fibers
d_fibers = [2 4 6 8 10 12 14];
dz = 5; % um
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
    for j=1:length(y_els)
        y_el = y_els(j);
        opts{end+1} = {"uNG"            2      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     y_el};
        opts{end+1} = {"TIME"           1      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     y_el};
    end
end
% unmyelinated fibers
d_fibers = [0.1 0.5 1 1.5 2];
z_shift = [0; 0.5];
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    for j=1:length(y_els)
        y_el = y_els(j);
        opts{end+1} = {"uNG"            2      4       200     0       nan         "ascent_sundt15"         d_fiber         0          z_shift     y_el};
        opts{end+1} = {"TIME"           1      4       200     0       nan         "ascent_sundt15"         d_fiber         0          z_shift     y_el};
    end
end
% =============================================================================================================================
% OVERLAP WINDOWS - EXTRAFASCICULAR (uNG, TIME)
y_els = -190-10-[2 50];
% myelinated fibers
d_fibers = [2 4 6 8 10 12 14];
dz = 5; % um
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
    for j=1:length(y_els)
        y_el = y_els(j);
        % opts_i{6} = 10 not used !!!
        opts{end+1} = {"uNG"         2      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     y_el};
        opts{end+1} = {"TIME"        1      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     y_el};
    end
end
% unmyelinated fibers
d_fibers = [0.1 0.5 1 1.5 2];
z_shift = [0; 0.5];
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    for j=1:length(y_els)
        y_el = y_els(j);
        % opts_i{6} = 10 not used !!!
        opts{end+1} = {"uNG"         2      4       200     0       10          "ascent_sundt15"         d_fiber         -190       z_shift     y_el};
        opts{end+1} = {"TIME"        1      4       200     0       10          "ascent_sundt15"         d_fiber         -190       z_shift     y_el};
    end
end
% =============================================================================================================================
% OVERLAP WINDOWS - EXTRANEURAL (cuff)
d_as_fascs = [20 500]; 
d_fasc_fib = 200;
% myelinated fibers
d_fibers = [2 4 6 8 10 12 14];
dz = 5; % um
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
    for j=1:length(d_as_fascs)
        d_as_fasc = d_as_fascs(j);
        opts{end+1} = {"Cortec_Cuff"        2      4       200     -800+d_as_fasc       nan         "ascent_mrg"         d_fiber         -1000+d_fasc_fib+d_as_fasc       z_shift     nan};
    end
end
% unmyelinated fibers
d_fibers = [0.1 0.5 1 1.5 2.0];
z_shift = [0; 0.5];
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    for j=1:length(d_as_fascs)
        d_as_fasc = d_as_fascs(j);
        opts{end+1} = {"Cortec_Cuff"        2      4       200     -800+d_as_fasc       nan         "ascent_sundt15"    d_fiber         -1000+d_fasc_fib+d_as_fasc       z_shift     nan};
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
                hm_params.electrodes.y = -hm_params.electrodes.h_shaft/2 + opts_i{11}*1e-6;
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
                hm_params.electrodes.y = opts_i{11}*1e-6;
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