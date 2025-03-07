%% PLOT FIGURES 3D, 4D, S3C, S4C
% the data must be simulated with the script run_simulations/sim2_SUAPs_along_fiber.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber");

%% PLOT OVERLAP WINDOW OF SUAP AMPLITUDES

d_fibers = [2:2:14];
d_fibers_unmyel = [0.1 0.5 1 1.5 2.0];

% uncomment one of the code blocks below to plot the desired figure
% =============================================================================================================================
% INTRAFASCICULAR TIME @ 2-50 um from fiber
elec = "TIME";
elec_mesh = 1;
y_els = -[2 50];
Vpp_lim_violins = [0 75];
dz = 5; % um
median_colors = [];
fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % INTRAFASCICULAR uNG @ 2-50 um from fiber
% elec = "uNG";
% elec_mesh = 2;
% y_els = -[2 50];
% Vpp_lim_violins = [0 75];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR TIME @ 2-50 um from fascicle
% elec = "TIME";
% elec_mesh = 1;
% y_els = -190-10-[2 50];
% Vpp_lim_violins = [0 15];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_mrg"    d_fiber      -190            z_shift     y_el};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_sundt15"    d_fiber     -190         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR uNG @ 2-50 um from fascicle
% elec = "uNG";
% elec_mesh = 2;
% y_els = -190-10-[2 50];
% Vpp_lim_violins = [0 7];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_mrg"    d_fiber      -190            z_shift     y_el};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10         "ascent_sundt15"    d_fiber     -190         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRANEURAL CUFF @ 20-500 um from fascicle
% elec = "Cortec_Cuff";
% elec_mesh = 2;
% Vpp_lim_violins = [0 18];
% % Vpp_lim_violins = [0 0.165]; % close-up
% d_as_fasc = [20 500];
% y_els = -800+d_as_fasc; % NOTE: is the y-position of the fascicle center!
% d_fasc_fib = 200;
% dz = 5; % um
% median_colors = colors_diams;
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     y_el        nan        "ascent_mrg"    d_fiber      y_el-200+d_fasc_fib         z_shift     nan};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     y_el        nan        "ascent_sundt15"    d_fiber      y_el-200+d_fasc_fib         z_shift_unmyel     nan};
% =============================================================================================================================

colors_diams = lines();
% colors_diams(3:4,:) = [];

% retrieve Vpp data - myelinated fibers
z_shift = cell(length(d_fibers), 1);
N1 = cell(length(d_fibers), length(y_els));
TN1 = cell(length(d_fibers), length(y_els));
P1 = cell(length(d_fibers), length(y_els));
TP1 = cell(length(d_fibers), length(y_els));
Vpp = cell(length(d_fibers), length(y_els));
Tpp = cell(length(d_fibers), length(y_els));
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift{i}, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);

    for j=1:length(y_els)
        y_el = y_els(j);
        opts = fn_opts(elec, elec_mesh, d_fiber, z_shift{i}, y_el);

        dir_model = fullfile(dir_experiments, DataHash(opts, 'MD5'));

        if strcmp(elec,'Cortec_Cuff')
            SUAP_t = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_t;
            SUAP_dt = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_dt;
            SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
            SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
            SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
        else
            load(fullfile(dir_model,"SUAP_1_1.mat"));
        end

        N1{i,j} = zeros(size(z_shift{i}));
        TN1{i,j} = zeros(size(z_shift{i}));
        P1{i,j} = zeros(size(z_shift{i}));
        TP1{i,j} = zeros(size(z_shift{i}));
        Vpp{i,j} = zeros(size(z_shift{i}));
        Tpp{i,j} = zeros(size(z_shift{i}));
        for n=1:length(z_shift{i})
            if strcmp(elec,'Cortec_Cuff')
                SUAP_n = (2*SUAP2{n}-SUAP1{n}-SUAP3{n})/2;
            else
                SUAP_n = SUAP{n};
            end

            % Vpp{i,j}(n) = max(SUAP_n)-min(SUAP_n);
    
            [N1{i,j}(n), TN1_idx] = min(SUAP_n);
            TN1{i,j}(n) = SUAP_t{n}(TN1_idx);
            [P1{i,j}(n), TP1_idx] = max(SUAP_n(SUAP_t{n}<TN1{i,j}(n)));
            TP1{i,j}(n) = SUAP_t{n}(TP1_idx);
            Vpp{i,j}(n) = P1{i,j}(n)-N1{i,j}(n);
            Tpp{i,j}(n) = TN1{i,j}(n)-TP1{i,j}(n);
        end
    end
end

% retrieve Vpp data - unmyelinated fibers
N1_unmyel = cell(length(d_fibers_unmyel), length(y_els));
TN1_unmyel = cell(length(d_fibers_unmyel), length(y_els));
P1_unmyel = cell(length(d_fibers_unmyel), length(y_els));
TP1_unmyel = cell(length(d_fibers_unmyel), length(y_els));
Vpp_unmyel = cell(length(d_fibers_unmyel), length(y_els));
Tpp_unmyel = cell(length(d_fibers_unmyel), length(y_els));
for i=1:length(d_fibers_unmyel)
    d_fiber = d_fibers_unmyel(i);
    z_shift_unmyel = [0; 0.5];

    for j=1:length(y_els)
        y_el = y_els(j);
        opts = fn_opts_unmyel(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el);

        dir_model = fullfile(dir_experiments, DataHash(opts, 'MD5'));

        if strcmp(elec,'Cortec_Cuff')
            SUAP_t = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_t;
            SUAP_dt = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_dt;
            SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
            SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
            SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
        else
            load(fullfile(dir_model,"SUAP_1_1.mat"));
        end

        N1_unmyel{i,j} = zeros(size(z_shift_unmyel));
        TN1_unmyel{i,j} = zeros(size(z_shift_unmyel));
        P1_unmyel{i,j} = zeros(size(z_shift_unmyel));
        TP1_unmyel{i,j} = zeros(size(z_shift_unmyel));
        Vpp_unmyel{i,j} = zeros(size(z_shift_unmyel));
        Tpp_unmyel{i,j} = zeros(size(z_shift_unmyel));
        for n=1:length(z_shift_unmyel)
            if strcmp(elec,'Cortec_Cuff')
                SUAP_n = (2*SUAP2{n}-SUAP1{n}-SUAP3{n})/2;
            else
                SUAP_n = SUAP{n};
            end

            % Vpp_unmyel{i,j}(n) = max(SUAP_n)-min(SUAP_n);

            [N1_unmyel{i,j}(n), TN1_idx] = min(SUAP_n);
            TN1_unmyel{i,j}(n) = SUAP_t{n}(TN1_idx);
            [P1_unmyel{i,j}(n), TP1_idx] = max(SUAP_n(SUAP_t{n}<TN1_unmyel{i,j}(n)));
            TP1_unmyel{i,j}(n) = SUAP_t{n}(TP1_idx);
            Vpp_unmyel{i,j}(n) = P1_unmyel{i,j}(n)-N1_unmyel{i,j}(n);
            Tpp_unmyel{i,j}(n) = TN1_unmyel{i,j}(n)-TP1_unmyel{i,j}(n);
        end
    end
end

% plot violinplots
for j=1:length(y_els)
    y_el = y_els(j);

    % select variables to plot
    data_myel = Vpp;
    data_unmyel = Vpp_unmyel;
    
    % data_myel = Tpp;
    % data_unmyel = Tpp_unmyel;
    
    pos_rel = 1:length(d_fibers);
    width_vs = 0.8;

    figure; 
    for i=1:length(d_fibers)
        if isempty(median_colors)
            median_color = [1 1 1];
        else
            median_color = median_colors(i,:);
        end
        vs(i) = Violin(data_myel(i,j),pos_rel(i),...
            'DataStyle', 'none',... % scatter, histogram
            'ShowMedian', true,...
            'MedianMarkerSize', 70, ...
            'MedianColor', median_color, ...
            'ShowWhiskers', true,...
            'QuartileStyle', 'shadow',...
            'Width', width_vs, ...
            'ViolinColor', {colors_diams(i,:)}, ...
            "HalfViolin",'left');
        h_lines(i) = plot(nan, nan, '-', 'Color', colors_diams(i,:), 'LineWidth', 1);
    end
    if iscell(data_unmyel)
        data_unmyel = arrayfun(@(x) x{1}(1), data_unmyel);
    end
    xlim([0 length(d_fibers)+1-width_vs]);
    for k=1:length(d_fibers_unmyel)
        h_lines_unmyel(k) = plot(xlim, data_unmyel(k,j)*[1 1], '-.', ...
            'Color', colors_diams(k,:), 'LineWidth', 1);
    end
    if ~isempty(Vpp_lim_violins)
        ylim(Vpp_lim_violins);
    end
    set(gca, 'FontSize', 15);
    set(gca, 'XTickLabel', [], 'XTick', []);
    ylabel('Vpp [uV]', 'FontSize', 15);
    legend([h_lines h_lines_unmyel], ["Myelinated "+d_fibers+" um" "Unmyelinated "+d_fibers_unmyel+" um"], 'Location', 'EastOutside');
    set(gcf, 'units', 'centimeters');
    pos = get(gcf,'Position');
    set(gcf, 'Position', [pos(1:2) 17 11]);
end