%% PLOT FIGURES 6D
% the data must be simulated with the script run_simulations/sim5_SUAPs_along_fiber_cuff_design.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber_cuff_design");


% parameters to retrieve data
R_cuff = [0.25 0.5 1];
L_cuff = [1 1.1 1.25 1.5 2 3 4 5 7.5 10 12.5 15 20];
d_fibers = [4 8];
dz = 5; % um
d_fibers_unmyel = [1 1.5];

% retrieve Vpp data - myelinated fibers
z_shift = cell(length(d_fibers), 1);
Vpp = cell(length(d_fibers),1);
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift{i}, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);

    Vpp{i} = zeros(length(R_cuff), length(L_cuff), length(z_shift{i}));
    for j=1:length(R_cuff)
        for k=1:length(L_cuff)
            opts = {"Cortec_Cuff"   2      4       200     0       nan        "ascent_mrg"         d_fiber      0      z_shift{i}     nan      R_cuff(j)   L_cuff(k)};
            dir_model = fullfile(dir_experiments, DataHash(opts, 'MD5'));
    
            load(fullfile(dir_model,"SUAP_1_1.mat"));
            for n=1:length(z_shift{i})
                Vpp{i}(j,k,n) = max(SUAP{n})-min(SUAP{n});
            end
        end
    end
end

% retrieve Vpp data - unmyelinated fibers
Vpp_unmyel = cell(length(d_fibers_unmyel), 1);
for i=1:length(d_fibers_unmyel)
    d_fiber_unmyel = d_fibers_unmyel(i);
    z_shift_unmyel = [0; 0.5];

    Vpp_unmyel{i} = zeros(length(R_cuff), length(L_cuff));
    for j=1:length(R_cuff)
        for k=1:length(L_cuff)
            opts = {"Cortec_Cuff"   2      4       200     0       nan        "ascent_sundt15"     d_fiber_unmyel      0      z_shift_unmyel     nan     R_cuff(j)      L_cuff(k)};
            dir_model = fullfile(dir_experiments, DataHash(opts, 'MD5'));
    
            load(fullfile(dir_model,"SUAP_1_1.mat"));
            Vpp_unmyel{i}(j,k) = max(SUAP{1})-min(SUAP{1});
        end
    end
end

R_cuff_plot = [0.25 0.5 1.0];
R_cuff_plot_idcs = find(any(R_cuff_plot == R_cuff(:),2));
line_colors_palette = lines();

line_style = ["-o" "-o"];
line_colors = line_colors_palette([1 2], :);
for j=1:length(d_fibers)
    fig = figure;
    set(fig, 'units', 'centimeters');
    fig.Position([3 4]) = [10 11];
    hold on;
    for i=1:length(R_cuff_plot)
        plot(L_cuff, median(Vpp{j}(R_cuff_plot_idcs(i),:,:),3), line_style(j), ...
            'Color', line_colors(j,:), 'MarkerFaceColor', line_colors(j,:) , 'LineWidth', 1);
    end
    ylim([0 60]);
    xlim([1.5 max(L_cuff)]);
    xlabel('Lcuff [mm]');
    ylabel('Vpp [uV]');
    set(gca, 'FontSize', 15);
end

line_style = ["-." "-."];
line_colors = line_colors_palette([2 5], :);
for j=1:length(d_fibers_unmyel)
    fig = figure;
    set(fig, 'units', 'centimeters');
    fig.Position([3 4]) = [10 11];
    hold on;
    for i=1:length(R_cuff_plot)
        plot(L_cuff, Vpp_unmyel{j}(R_cuff_plot_idcs(i),:), line_style(j), ...
            'Color', line_colors(j,:), 'MarkerFaceColor', line_colors(j,:) , 'LineWidth', 1);
    end
    ylim([0 0.6]);
    xlim([1.5 max(L_cuff)]);
    xlabel('Lcuff [mm]');
    ylabel('Vpp [uV]');
    set(gca, 'FontSize', 15);
end
