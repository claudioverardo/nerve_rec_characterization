%% PLOT FIGURES 5E-F
% the data must be simulated with the script run_simulations/sim3_SUAPs_along_fiber_perineurium.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber_perineurium");


T_peri = [0.25 0.5 0.75 1.0 1.25] * FascicleTopography.perineurium_thickness("pelot20_human", 400e-6)*1e+6;
d_fibers = [4 8];
d_fibers_unmyel = [0.5 1.0 1.5];

% uncomment one of the code blocks below to plot the desired figure
% =============================================================================================================================
% EXTRAFASCICULAR TIME @ 2 um from fascicle
elec = "TIME";
elec_mesh = 1;
y_el = -190-10-[2];
Vpp_lim_violins = [0 30];
dz = 5; % um
median_colors = [];
fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      -190            z_shift     y_el     T_peri};
fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
    {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber     -190         z_shift     y_el    T_peri};
% =============================================================================================================================
% % EXTRAFASCICULAR uNG @ 2 um from fascicle
% elec = "uNG";
% elec_mesh = 2;
% y_el = -190-10-[2];
% Vpp_lim_violins = [0 8];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      -190            z_shift     y_el     T_peri};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber     -190         z_shift     y_el    T_peri};
% =============================================================================================================================
% % EXTRAFASCICULAR TIME @ 50 um from fascicle
% elec = "TIME";
% elec_mesh = 1;
% y_el = -190-10-[50];
% Vpp_lim_violins = [0 30];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      -190            z_shift     y_el     T_peri};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber     -190         z_shift     y_el    T_peri};
% =============================================================================================================================
% % EXTRAFASCICULAR uNG @ 50 um from fascicle
% elec = "uNG";
% elec_mesh = 2;
% y_el = -190-10-[50];
% Vpp_lim_violins = [0 8];
% dz = 5; % um
% median_colors = [];
% fn_opts = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_mrg"    d_fiber      -190            z_shift     y_el     T_peri};
% fn_opts_unmyel = @(elec, elec_mesh, d_fiber, z_shift, y_el, T_peri) ...
%     {elec   elec_mesh      4       200     0           nan        "ascent_sundt15"    d_fiber     -190         z_shift     y_el    T_peri};
% =============================================================================================================================

% retrieve Vpp data - myelinated fibers
z_shift = cell(length(d_fibers), 1);
N1 = cell(length(d_fibers), length(T_peri));
TN1 = cell(length(d_fibers), length(T_peri));
P1 = cell(length(d_fibers), length(T_peri));
TP1 = cell(length(d_fibers), length(T_peri));
Vpp = cell(length(d_fibers), length(T_peri));
Tpp = cell(length(d_fibers), length(T_peri));
for i=1:length(d_fibers)
    d_fiber = d_fibers(i);
    [z_shift{i}, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);

    for j=1:length(T_peri)
        opts = fn_opts(elec, elec_mesh, d_fiber, z_shift{i}, y_el, T_peri(j));
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
                SUAP_n = (2*SUAP1{n}-SUAP2{n}-SUAP3{n})/2;
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
N1_unmyel = cell(length(d_fibers_unmyel), length(T_peri));
TN1_unmyel = cell(length(d_fibers_unmyel), length(T_peri));
P1_unmyel = cell(length(d_fibers_unmyel), length(T_peri));
TP1_unmyel = cell(length(d_fibers_unmyel), length(T_peri));
Vpp_unmyel = cell(length(d_fibers_unmyel), length(T_peri));
Tpp_unmyel = cell(length(d_fibers_unmyel), length(T_peri));
for i=1:length(d_fibers_unmyel)
    d_fiber_unmyel = d_fibers_unmyel(i);
    z_shift_unmyel = [0; 0.5];

    for j=1:length(T_peri)
        opts = fn_opts_unmyel(elec, elec_mesh, d_fiber_unmyel, z_shift_unmyel, y_el, T_peri(j));
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
                SUAP_n = (2*SUAP1{n}-SUAP2{n}-SUAP3{n})/2;
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

% data_myel = Tpp;
% data_unmyel = Tpp_unmyel;
    
colors_diams = lines();
colors_diams(3:4,:) = [];
pos_rel = 1:length(T_peri);

fig = figure; 
set(fig, 'units', 'centimeters');
fig.Position([3]) = [11];
for j=1:length(T_peri)

    half_violin = ["left" "right"];
    width_vs = 0.4;
    shift_vs = 0.05;

    for i=1:length(d_fibers)
        if isempty(median_colors)
            median_color = [1 1 1];
        else
            median_color = median_colors(i,:);
        end
        vs(i) = Violin(data_myel(i,j),pos_rel(j)+(-1)^(i)*shift_vs,...
            'DataStyle', 'none',... % scatter, histogram
            'ShowMedian', true,...
            'MedianMarkerSize', 70, ...
            'MedianColor', median_color, ...
            'ShowWhiskers', true,...
            'QuartileStyle', 'shadow',...
            'Width', width_vs, ...
            'ViolinColor', {colors_diams(i,:)}, ...
            "HalfViolin",half_violin(mod(i+1,2)+1)); 
    end
end
if iscell(data_unmyel)
    data_unmyel = arrayfun(@(x) x{1}(1), data_unmyel);
end
for k=1:length(d_fibers_unmyel)
    h_lines_unmyel(k) = plot(pos_rel, data_unmyel(k,:), '-.', ...
        'Color', colors_diams(k,:), 'LineWidth', 1);
end
xlim([pos_rel(1)-width_vs-shift_vs-0.1 pos_rel(end)+width_vs+shift_vs+0.1]);
if ~isempty(Vpp_lim_violins)
    ylim(Vpp_lim_violins);
end
xlabel('T peri');
ylabel('Vpp [uV]');
set(gca, 'XTick', pos_rel, 'XTickLabel', ["0.25x" "0.5x" "0.75x" "1.0x" "1.25x"])
set(gca, 'FontSize', 15);
