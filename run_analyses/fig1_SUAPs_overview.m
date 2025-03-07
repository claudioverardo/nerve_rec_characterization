%% PLOT FIGURES 1G
% the data must be simulated with the script run_simulations/sim1_SUAPs.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis");

%% PLOT SUAPs (superimposed) - large-myelinated fibers
% @ half internode / 20 um distance from fiber

% uncomment the code blocks below to plot desired subset of figures
opts = {};
% ======================================================================================
% intrafascicular TIME - myelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     20          0.5};
% ======================================================================================
% intrafascicular uNG - myelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     20          0.5};
% ======================================================================================
% extrafascicular TIME - myelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       10          "ascent_mrg"         [1:0.5:13]     -200+20-10  0.5};
% ======================================================================================
% extrafascicular uNG - myelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       10          "ascent_mrg"         [1:0.5:13]     -200+20-10  0.5};
% ======================================================================================
% extraneural cuff - myelinated fibers
opts{end+1} = {"Cortec_Cuff"   2      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     0           0.5};
% ======================================================================================

diameters_lim = [6 10];
xlim_vec = [-0.8 1.0];
ylim_vec = [-10 8];

for n=1:length(opts)
    dir_model = fullfile(dir_experiments, DataHash(opts{n}, 'MD5'));

    if strcmp(opts{n}{1},'Cortec_Cuff')
        SUAP_t = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_t;
        SUAP_t0 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_t0;
        SUAP_CV = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_CV;
        SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
        SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
        SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
    else
        load(fullfile(dir_model,"SUAP_1_1.mat"));
    end

    diameters = opts{n}{8};
    diameters_sel = diameters(diameters>=diameters_lim(1) & diameters<=diameters_lim(2));
    CLim = [diameters_lim(1) diameters_lim(2)+1];
    colors = data2color(diameters_sel, CLim, copper(256));
    figure;
    for i=1:length(diameters_sel)
        idx = find(diameters_sel(i) == diameters);
        if strcmp(opts{n}{1},'Cortec_Cuff')
            SUAP_idx = (2*SUAP2{idx}-SUAP1{idx}-SUAP3{idx})/2;
        else
            SUAP_idx = SUAP{idx};
        end
        % RECORDING ELECTRODE ASSUMED AT z=0 !!!!!
        SUAP_t0_at_rec = SUAP_t0(idx) + (opts{n}{3}/2*1e-2)/SUAP_CV(idx)*1e+3;
        plot(SUAP_t{idx}-SUAP_t0_at_rec, SUAP_idx, 'Color', colors(i,:), 'LineWidth', 1);
        hold on;
    end
    xlim(xlim_vec);
    ylim(ylim_vec);
    h_cb = colorbar;
    h_cb.Limits = diameters_lim;
    clim(CLim);
    colormap(copper(256));
    set(gca, 'color', 'none', 'FontSize', 15);
    set(gca, 'FontName', 'Arial');
    set(h_cb, 'FontName', 'Arial');
    set(gcf, 'units', 'centimeters');
    pos = get(gcf,'Position');
    set(gcf, 'Position', [pos(1:2) 9 8]);
end

%% PLOT SUAPs (superimposed) - small-myelinated fibers
% @ half internode / 20 um distance from fiber

% uncomment the code blocks below to plot desired subset of figures
opts = {};
% ======================================================================================
% intrafascicular TIME - myelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     20          0.5};
% ======================================================================================
% intrafascicular uNG - myelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     20          0.5};
% ======================================================================================
% extrafascicular TIME - myelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       10          "ascent_mrg"         [1:0.5:13]     -200+20-10  0.5};
% ======================================================================================
% extrafascicular uNG - myelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       10          "ascent_mrg"         [1:0.5:13]     -200+20-10  0.5};
% ======================================================================================
% extraneural cuff - myelinated fibers
opts{end+1} = {"Cortec_Cuff"   2      4       200     0       nan         "ascent_mrg"         [1:0.5:13]     0           0.5};
% ======================================================================================

diameters_lim = [2 5];
xlim_vec = [-0.8 1.0];
ylim_vec = [-10 8];

for n=1:length(opts)
    dir_model = fullfile(dir_experiments, DataHash(opts{n}, 'MD5'));

    if strcmp(opts{n}{1},'Cortec_Cuff')
        SUAP_t = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_t;
        SUAP_t0 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_t0;
        SUAP_CV = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_CV;
        SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
        SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
        SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
    else
        load(fullfile(dir_model,"SUAP_1_1.mat"));
    end

    diameters = opts{n}{8};
    diameters_sel = diameters(diameters>=diameters_lim(1) & diameters<=diameters_lim(2));
    CLim = [diameters_lim(1) diameters_lim(2)+1];
    colors = data2color(diameters_sel, CLim, copper(256));
    figure;
    for i=1:length(diameters_sel)
        idx = find(diameters_sel(i) == diameters);
        if strcmp(opts{n}{1},'Cortec_Cuff')
            SUAP_idx = (2*SUAP2{idx}-SUAP1{idx}-SUAP3{idx})/2;
        else
            SUAP_idx = SUAP{idx};
        end
        % RECORDING ELECTRODE ASSUMED AT z=0 !!!!!
        SUAP_t0_at_rec = SUAP_t0(idx) + (opts{n}{3}/2*1e-2)/SUAP_CV(idx)*1e+3;
        plot(SUAP_t{idx}-SUAP_t0_at_rec, SUAP_idx, 'Color', colors(i,:), 'LineWidth', 1);
        hold on;
    end
    xlim(xlim_vec);
    ylim(ylim_vec);
    h_cb = colorbar;
    h_cb.Limits = diameters_lim;
    clim(CLim);
    colormap(copper(256));
    set(gca, 'color', 'none', 'FontSize', 15);
    set(gca, 'FontName', 'Arial');
    set(h_cb, 'FontName', 'Arial');
    set(gcf, 'units', 'centimeters');
    pos = get(gcf,'Position');
    set(gcf, 'Position', [pos(1:2) 9 8]);
end

%% PLOT SUAPs (superimposed) - unmyelinated fibers
% @ half internode / 20 um distance from fiber

% uncomment the code blocks below to plot desired subset of figures
opts = {};
% ======================================================================================
% intrafascicular TIME - unmyelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       nan         "ascent_sundt15"     [0.1:0.1:1.5]  20          0.5};
% ======================================================================================
% intrafascicular uNG - unmyelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       nan         "ascent_sundt15"     [0.1:0.1:1.5]  20          0.5};
% ======================================================================================
% extrafascicular TIME - unmyelinated fibers
opts{end+1} = {"TIME"          1      4       200     0       10          "ascent_sundt15"     [0.1:0.1:1.5]  -200+20-10  0.5};
% ======================================================================================
% extrafascicular uNG - unmyelinated fibers
opts{end+1} = {"uNG"           2      4       200     0       10          "ascent_sundt15"     [0.1:0.1:1.5]  -200+20-10  0.5};
% ======================================================================================
% extraneural cuff - unmyelinated fibers
opts{end+1} = {"Cortec_Cuff"   2      4       200     0       nan         "ascent_sundt15"     [0.1:0.1:1.5]  0           0.5};
% ======================================================================================

diameters_lim = [0.1 1.0];
xlim_vec = [-1.3 2.3];
ylim_vec = [-10 8];

for n=1:length(opts)
    dir_model = fullfile(dir_experiments, DataHash(opts{n}, 'MD5'));
    
    if strcmp(opts{n}{1},'Cortec_Cuff')
        SUAP_t = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_t;
        SUAP_t0 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_t0;
        SUAP_CV = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP_CV;
        SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
        SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
        SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
    else
        load(fullfile(dir_model,"SUAP_1_1.mat"));
    end

    diameters = opts{n}{8};
    diameters_sel = diameters(diameters>=diameters_lim(1) & diameters<=diameters_lim(2));
    CLim = [diameters_lim(1) diameters_lim(2)+0.2];
    colors = data2color(diameters_sel, CLim, copper(256));
    figure;
    for i=1:length(diameters_sel)
        idx = find(diameters_sel(i) == diameters);
        if strcmp(opts{n}{1},'Cortec_Cuff')
            SUAP_idx = (2*SUAP2{idx}-SUAP1{idx}-SUAP3{idx})/2;
        else
            SUAP_idx = SUAP{idx};
        end
        % RECORDING ELECTRODE ASSUMED AT z=0 !!!!!
        SUAP_t0_at_rec = SUAP_t0(idx) + (opts{n}{3}/2*1e-2)/SUAP_CV(idx)*1e+3;
        plot(SUAP_t{idx}-SUAP_t0_at_rec, SUAP_idx, 'Color', colors(i,:), 'LineWidth', 1);
        hold on;
    end
    xlim(xlim_vec);
    ylim(ylim_vec);
    h_cb = colorbar;
    h_cb.Limits = diameters_lim;
    clim(CLim);
    colormap(copper(256));
    set(gca, 'color', 'none', 'FontSize', 15);
    set(gca, 'FontName', 'Arial');
    set(h_cb, 'FontName', 'Arial');
    set(gcf, 'units', 'centimeters');
    pos = get(gcf,'Position');
    set(gcf, 'Position', [pos(1:2) 9 8]);
end

