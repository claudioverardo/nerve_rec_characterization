%% PLOT FIGURES 3A-C, 4A-C, S3A-B, S4A-B
% the data must be simulated with the script run_simulations/sim2_SUAPs_along_fiber.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber");

%% PLOT SUAP WAVEFORMS ALONG FIBERS @ FIXED ELECTRODE-FIBER DISTANCE
% - SUAPs are plotted for selected diameters of myelinated and unmyelinated fibers
% - SUAPs are plotted at Ranvier nodes, 1/4 of internode, and 1/2 of internode

d_fibers = [4 8];
z_shift_myel = [1.0 0.75 0.5];
d_fibers_unmyel = [1.0 1.5];
z_shift_unmyel = [0; 0.5];

% uncomment one of the code blocks below to plot the desired figure
% =============================================================================================================================
% INTRAFASCICULAR TIME @ 2 um from fiber
elec = "TIME";
elec_mesh = 1;
y_els = -[2];
Vpp_lim_SUAPs_myel = [-61 30];
t_lim_SUAPs_myel = [-0.2 0.4];
dz = 5; % um
fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
Vpp_lim_SUAPs_unmyel = [-61 30];
t_lim_SUAPs_unmyel = [-1.1 2.1];
fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % INTRAFASCICULAR uNG @ 2 um from fiber
% elec = "uNG";
% elec_mesh = 2;
% y_els = -[2];
% Vpp_lim_SUAPs_myel = [-52 22];
% t_lim_SUAPs_myel = [-0.2 0.4];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
% Vpp_lim_SUAPs_unmyel = [-52 22];
% t_lim_SUAPs_unmyel = [-1.1 2.1];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR TIME @ 2 um from fascicle
% elec = "TIME";
% elec_mesh = 1;
% y_els = -190-10-[2];
% Vpp_lim_SUAPs_myel = [-9 5];
% t_lim_SUAPs_myel = [-0.2 0.4];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_mrg"    d_fiber      -190            z_shift     y_el};
% Vpp_lim_SUAPs_unmyel = [-9 5];
% t_lim_SUAPs_unmyel = [-1.1 2.1];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_sundt15"    d_fiber     -190         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR uNG @ 2 um from fascicle
% elec = "uNG";
% elec_mesh = 2;
% y_els = -190-10-[2];
% Vpp_lim_SUAPs_myel = [-2 1.2];
% t_lim_SUAPs_myel = [-0.2 0.4];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_mrg"    d_fiber      -190            z_shift     y_el};
% Vpp_lim_SUAPs_unmyel = [-2 1.2];
% t_lim_SUAPs_unmyel  = [-1.1 2.1];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_sundt15"    d_fiber     -190         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRANEURAL CUFF @ 20 um from fascicle
% elec = "Cortec_Cuff";
% elec_mesh = 2;
% d_as_fasc = [20];
% y_els = -800+d_as_fasc; % NOTE: is the y-position of the fascicle center!
% d_fasc_fib = 200;
% Vpp_lim_SUAPs_myel = [-4.5 3.5];
% t_lim_SUAPs_myel = [-0.4 0.6];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     y_el        nan        "ascent_mrg"    d_fiber      y_el-200+d_fasc_fib         z_shift     nan};
% Vpp_lim_SUAPs_unmyel = [-4.5 3.5];
% t_lim_SUAPs_unmyel = [-17 17];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     y_el        nan        "ascent_sundt15"    d_fiber      y_el-200+d_fasc_fib         z_shift_unmyel     nan};
% =============================================================================================================================
% % INTRAFASCICULAR TIME @ 50 um from fiber
% elec = "TIME";
% elec_mesh = 1;
% y_els = -[50];
% Vpp_lim_SUAPs_myel = [-21 12];
% t_lim_SUAPs_myel = [-0.2 0.4];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
% Vpp_lim_SUAPs_unmyel = [-21 12];
% t_lim_SUAPs_unmyel = [-1.1 2.1];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % INTRAFASCICULAR uNG @ 50 um from fiber
% elec = "uNG";
% elec_mesh = 2;
% y_els = -[50];
% Vpp_lim_SUAPs_myel = [-13 8];
% t_lim_SUAPs_myel = [-0.2 0.4];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
% Vpp_lim_SUAPs_unmyel = [-13 8];
% t_lim_SUAPs_unmyel = [-1.1 2.1];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR TIME @ 50 um from fascicle
% elec = "TIME";
% elec_mesh = 1;
% y_els = -190-10-[50];
% Vpp_lim_SUAPs_myel = [-4 2.5];
% t_lim_SUAPs_myel = [-0.2 0.4];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_mrg"    d_fiber      -190            z_shift     y_el};
% Vpp_lim_SUAPs_unmyel = [-4 2.5];
% t_lim_SUAPs_unmyel = [-1.1 2.1];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_sundt15"    d_fiber     -190         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR uNG @ 50 um from fascicle
% elec = "uNG";
% elec_mesh = 2;
% y_els = -190-10-[50];
% Vpp_lim_SUAPs_myel = [-2 1.2];
% t_lim_SUAPs_myel = [-0.2 0.4];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_mrg"    d_fiber      -190            z_shift     y_el};
% Vpp_lim_SUAPs_unmyel = [-2 1.2];
% t_lim_SUAPs_unmyel  = [-1.1 2.1];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_sundt15"    d_fiber     -190         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRANEURAL CUFF @ 500 um from fascicle
% elec = "Cortec_Cuff";
% elec_mesh = 2;
% d_as_fasc = [500];
% y_els = -800+d_as_fasc; % NOTE: is the y-position of the fascicle center!
% d_fasc_fib = 200;
% Vpp_lim_SUAPs_myel = [-4.5 3.5];
% t_lim_SUAPs_myel = [-0.4 0.6];
% dz = 5; % um
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     y_el        nan        "ascent_mrg"    d_fiber      y_el-200+d_fasc_fib         z_shift     nan};
% Vpp_lim_SUAPs_unmyel = [-4.5 3.5];
% t_lim_SUAPs_unmyel = [-17 17];
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     y_el        nan        "ascent_sundt15"    d_fiber      y_el-200+d_fasc_fib         z_shift_unmyel     nan};
% =============================================================================================================================

colors_lines = lines();
colors_lines(3:4,:) = [];

% plot myelinated SUAPS
for j=1:length(y_els)
    y_el = y_els(j);

    for i=1:length(d_fibers)
        figure; 

        d_fiber = d_fibers(i);
        [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
        opts = fn_opts(elec, elec_mesh, d_fiber, z_shift, y_el);
        dir_model = fullfile(dir_experiments, DataHash(opts, 'MD5'));
        
        for n=1:length(z_shift_myel)
            [~, idx] = min(abs(z_shift-z_shift_myel(n)));
            if strcmp(elec,'Cortec_Cuff')
                SUAP_t = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_t;
                SUAP_t0 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_t0;
                SUAP_CV = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_CV;
                SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
                SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
                SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
                SUAP_idx = (2*SUAP2{idx}-SUAP1{idx}-SUAP3{idx})/2;
            else
                load(fullfile(dir_model,"SUAP_1_1.mat"));
                SUAP_idx = SUAP{idx};
            end
    
            % RECORDING ELECTRODE ASSUMED AT z=0 !!!!!
            SUAP_t0_at_rec = SUAP_t0(idx) + (opts{3}/2*1e-2)/SUAP_CV(idx)*1e+3;
            plot(SUAP_t{idx}-SUAP_t0_at_rec, SUAP_idx, 'Color', colors_lines(n,:), 'LineWidth', 1);
            hold on;
        end
        box off;
        ylim(Vpp_lim_SUAPs_myel);
        xlim(t_lim_SUAPs_myel);
        ylabel('SUAP [uV]', 'FontSize', 15);
        xlabel('Time [ms]', 'FontSize', 15);
        set(gca, 'FontSize', 15);
        set(gcf, 'units', 'centimeters');
        pos = get(gcf,'Position');
        set(gcf, 'Position', [pos(1:2) 8 11]);
    end
end

% plot unmyelinated SUAPS
for j=1:length(y_els)
    y_el = y_els(j);

    for k=1:length(d_fibers_unmyel)
        figure; 

        d_fiber = d_fibers_unmyel(k);
        opts = fn_opts_unmyel(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el);
        dir_model = fullfile(dir_experiments, DataHash(opts, 'MD5'));

        idx = 1;
        if strcmp(elec,'Cortec_Cuff')
            SUAP_t = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_t;
            SUAP_dt = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_dt;
            SUAP_t0 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_t0;
            SUAP_CV = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_CV;
            SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
            SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
            SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
            SUAP_idx = (2*SUAP2{idx}-SUAP1{idx}-SUAP3{idx})/2;
        else
            load(fullfile(dir_model,"SUAP_1_1.mat"));
            SUAP_idx = SUAP{idx};
        end

        % RECORDING ELECTRODE ASSUMED AT z=0 !!!!!
        SUAP_t0_at_rec = SUAP_t0(idx) + (opts{3}/2*1e-2)/SUAP_CV(idx)*1e+3;
        plot(SUAP_t{idx}-SUAP_t0_at_rec, SUAP_idx, '-', 'Color', colors_lines(1,:), 'LineWidth', 1);
        hold on;

        box off;
        ylim(Vpp_lim_SUAPs_unmyel);
        xlim(t_lim_SUAPs_unmyel);
        ylabel('SUAP [uV]', 'FontSize', 15);
        xlabel('Time [ms]', 'FontSize', 15);
        set(gca, 'FontSize', 15);
        set(gcf, 'units', 'centimeters');
        pos = get(gcf,'Position');
        set(gcf, 'Position', [pos(1:2) 8 11]);
    end

end

