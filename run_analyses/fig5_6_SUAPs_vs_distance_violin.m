%% PLOT FIGURES 5A-D, 6A
% the data must be simulated with the script run_simulations/sim2_SUAPs_along_fiber.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber");


d_fibers = [4 8];
d_fibers_unmyel = [0.5 1 1.5];
colors_diams = lines();
colors_diams = colors_diams([1 2 1 2 5],:);

% uncomment one of the code blocks below to plot the desired figure
% =============================================================================================================================
% INTRAFASCICULAR TIME
elec = "TIME";
elec_mesh = 1;
y_els = -[0:10:100];
y_els_unmyel = -[0:2:100];
data_ylim = [0 60];
data_xlim = [-5 105];
dz = 5; % um
median_colors = [];
fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % INTRAFASCICULAR uNG
% elec = "uNG";
% elec_mesh = 2;
% y_els = -[0:10:100];
% y_els_unmyel = -[0:2:100];
% data_ylim = [0 60];
% data_xlim = [-5 105];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      0        z_shift     y_el};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber      0         z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR TIME
% elec = "TIME";
% elec_mesh = 1;
% y_els = -190-10-[0:10:100];
% y_els_unmyel = -190-10-[0:2:100];
% data_ylim = [0 15];
% data_xlim = [195 305];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10        "ascent_mrg"    d_fiber      -190        z_shift           y_el};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10        "ascent_sundt15"    d_fiber   -190       z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRAFASCICULAR uNG
% elec = "uNG";
% elec_mesh = 2;
% y_els = -190-10-[0:10:100];
% y_els_unmyel = -190-10-[0:2:100];
% data_ylim = [0 4];
% data_xlim = [195 305];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     0           10        "ascent_mrg"    d_fiber      -190        z_shift           y_el};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     0           10        "ascent_sundt15"    d_fiber   -190       z_shift_unmyel     y_el};
% =============================================================================================================================
% % EXTRANEURAL CUFF
% elec = "Cortec_Cuff";
% elec_mesh = 2;
% y_els = -[0:20:200 300:100:800];
% y_els_unmyel = -[0:20:200 300:100:800];
% d_fasc_fib = 200;
% data_ylim = [0 6.5];
% % data_ylim = [0 0.1]; % close-up
% data_xlim = [-50 850];
% dz = 5; % um
% median_colors = colors_diams;
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el) ...
%     {elec   elec_mesh      4       200     -800-y_el   nan        "ascent_mrg"    d_fiber      -1000+d_fasc_fib-y_el   z_shift     nan};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el) ...
%     {elec   elec_mesh      4       200     -800-y_el   nan        "ascent_sundt15"    d_fiber      -1000+d_fasc_fib-y_el   z_shift_unmyel     nan};
% =============================================================================================================================

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
            load(fullfile(dir_model,"SUAP_1_1.mat"), 'SUAP_t', 'SUAP_dt');
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
N1_unmyel = cell(length(d_fibers_unmyel), length(y_els_unmyel));
TN1_unmyel = cell(length(d_fibers_unmyel), length(y_els_unmyel));
P1_unmyel = cell(length(d_fibers_unmyel), length(y_els_unmyel));
TP1_unmyel = cell(length(d_fibers_unmyel), length(y_els_unmyel));
Vpp_unmyel = cell(length(d_fibers_unmyel), length(y_els_unmyel));
Tpp_unmyel = cell(length(d_fibers_unmyel), length(y_els_unmyel));
for i=1:length(d_fibers_unmyel)
    d_fiber = d_fibers_unmyel(i);
    z_shift_unmyel = [0; 0.5];

    for j=1:length(y_els_unmyel)
        y_el = y_els_unmyel(j);
        opts = fn_opts_unmyel(elec, elec_mesh, d_fiber, z_shift_unmyel, y_el);
        dir_model = fullfile(dir_experiments, DataHash(opts, 'MD5'));

        if strcmp(elec,'Cortec_Cuff')
            load(fullfile(dir_model,"SUAP_1_1.mat"), 'SUAP_t', 'SUAP_dt');
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

% plot violins
data_myel = Vpp;
data_unmyel = Vpp_unmyel;
ylabel_str = 'Vpp across fiber [uV]';

% data_myel = Tpp;
% data_unmyel = Tpp_unmyel;
% ylabel_str = 'Tpp across fiber [ms]';

width_vs = 3;
pos_rel = -(length(d_fibers)-1)/2:(length(d_fibers)-1)/2;
half_violin = ["left" "right"];

figure; 
idx = 0;
for j=1:length(y_els)
    for i=1:length(d_fibers)
        idx = idx + 1;
        if isempty(median_colors)
            median_color = [1 1 1];
        else
            median_color = colors_diams(i,:);
        end
        vs(idx) = Violin(data_myel(i,j),-y_els(j)+pos_rel(i),...
            'DataStyle', 'none',... % scatter, histogram
            'ShowMedian', true,...
            'ShowWhiskers', true,...
            'MedianColor',  median_color,...
            'QuartileStyle', 'shadow',...
            'Width', width_vs, ...
            'ViolinColor', {colors_diams(i,:)}, ...
            "HalfViolin",half_violin(mod(i+1,2)+1));
    end
end
hold on;
if iscell(data_unmyel)
    data_unmyel = arrayfun(@(x) x{1}(1), data_unmyel);
end
for k=1:length(d_fibers_unmyel)
    h_lines(k) = plot(-y_els_unmyel, data_unmyel(k,:), '-.', 'Color', colors_diams(i+k,:), 'LineWidth', 1);
end

legend([vs(1:length(d_fibers)).ViolinPlotQ h_lines], ["Myelinated "+d_fibers+" um" "Unmyelinated "+d_fibers_unmyel+" um"]);
xlabel('Distance from fiber [um]', 'FontSize', 15);
ylabel(ylabel_str, 'FontSize', 15);
set(gca, 'FontSize', 15);
if ~isempty(data_xlim)
    xlim(data_xlim);
end
if ~isempty(data_ylim)
    ylim(data_ylim);
end
set(gcf, 'units', 'centimeters');
pos = get(gcf,'Position');
set(gcf, 'Position', [pos(1:2) 18 9]);
