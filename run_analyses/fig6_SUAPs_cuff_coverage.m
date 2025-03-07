%% PLOT FIGURES 6C
% the data must be simulated with the script run_simulations/sim4_SUAPs_along_fiber_cuff_coverage.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber_cuff_coverage");


% parameters to retrieve data
elec = "Cortec_Cuff";
elec_mesh = 2;
Vpp_lim = [0 11];
% Vpp_lim = [0 0.3]; % close-up
R_nerve = 1e-3; % m
H_AS = [0.02 0.1 0.2 0.5 1.0]*2*pi*R_nerve*1e+3;
W_AS = [1.0  1.0 1.0 1.0 1.0];
d_as_fasc = 20;
d_fasc_fib = 200;
d_fibers = [4 8];
d_fibers_unmyel = [0.5 1.0 1.5];

% retrieve Vpp data - myelinated fibers
N1 = zeros(length(d_fibers), length(H_AS));
TN1 = zeros(length(d_fibers), length(H_AS));
P1 = zeros(length(d_fibers), length(H_AS));
TP1 = zeros(length(d_fibers), length(H_AS));
Vpp = zeros(length(d_fibers), length(H_AS));
Tpp = zeros(length(d_fibers), length(H_AS));
for i=1:length(H_AS)
    opts = {elec   elec_mesh      4       200     -800+d_as_fasc   nan        "ascent_mrg"       [2 4 6 8 10 12 14]     -1000+d_fasc_fib+d_as_fasc   0   nan   H_AS(i)   W_AS(i)};
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

    for n=1:length(d_fibers)
        idx = find(opts{8} == d_fibers(n));
        if strcmp(elec,'Cortec_Cuff')
            SUAP_n = (2*SUAP2{idx}-SUAP1{idx}-SUAP3{idx})/2;
        else
            SUAP_n = SUAP{idx};
        end

        % Vpp(n,i) = max(SUAP_n)-min(SUAP_n);

        [N1(n,i), TN1_idx] = min(SUAP_n);
        TN1(n,i) = SUAP_t{n}(TN1_idx);
        [P1(n,i), TP1_idx] = max(SUAP_n(SUAP_t{n}<TN1(n,i)));
        TP1(n,i) = SUAP_t{n}(TP1_idx);
        Vpp(n,i) = P1(n,i)-N1(n,i);
        Tpp(n,i) = TN1(n,i)-TP1(n,i);
    end
end

% retrieve Vpp data - unmyelinated fibers
N1_unmyel = zeros(length(d_fibers_unmyel), length(H_AS));
TN1_unmyel = zeros(length(d_fibers_unmyel), length(H_AS));
P1_unmyel = zeros(length(d_fibers_unmyel), length(H_AS));
TP1_unmyel = zeros(length(d_fibers_unmyel), length(H_AS));
Vpp_unmyel = zeros(length(d_fibers_unmyel), length(H_AS));
Tpp_unmyel = zeros(length(d_fibers_unmyel), length(H_AS));
for i=1:length(H_AS)
    opts = {elec   elec_mesh      4       200     -800+d_as_fasc   nan        "ascent_sundt15"    d_fibers_unmyel      -1000+d_fasc_fib+d_as_fasc  0   nan   H_AS(i)   W_AS(i)};
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

    for n=1:length(d_fibers_unmyel)
        idx = find(opts{8} == d_fibers_unmyel(n));
        if strcmp(elec,'Cortec_Cuff')
            SUAP_n = (2*SUAP2{idx}-SUAP1{idx}-SUAP3{idx})/2;
        else
            SUAP_n = SUAP{idx};
        end

        % Vpp_unmyel(n,i) = max(SUAP_n)-min(SUAP_n);

        [N1_unmyel(n,i), TN1_idx] = min(SUAP_n);
        TN1_unmyel(n,i) = SUAP_t{n}(TN1_idx);
        [P1_unmyel(n,i), TP1_idx] = max(SUAP_n(SUAP_t{n}<TN1_unmyel(n,i)));
        TP1_unmyel(n,i) = SUAP_t{n}(TP1_idx);
        Vpp_unmyel(n,i) = P1_unmyel(n,i)-N1_unmyel(n,i);
        Tpp_unmyel(n,i) = TN1_unmyel(n,i)-TP1_unmyel(n,i);
    end
end

% plot
line_colors_palette = lines();
xticks_vec = 1:length(H_AS);

fig = figure; 
set(fig, 'units', 'centimeters');
fig.Position([3 4]) = [11 10];
hold on;
line_style = ["-o" "-o"];
line_colors = line_colors_palette([1 2], :);
for i=1:length(d_fibers)
    plot(xticks_vec, Vpp(i,:), line_style(i), ...
        'Color', line_colors(i,:), 'MarkerFaceColor', line_colors(i,:) , 'LineWidth', 1);
end
line_style = ["-." "-." "-."];
line_colors = line_colors_palette([1 2 5], :);
for j=1:length(d_fibers_unmyel)
    plot(xticks_vec, Vpp_unmyel(j,:), line_style(j), ...
        'Color', line_colors(j,:), 'MarkerFaceColor', line_colors(j,:) , 'LineWidth', 1);
end
xlim([xticks_vec(1)-0.4 xticks_vec(end)+0.4]);
if ~isempty(Vpp_lim)
    ylim(Vpp_lim);
end
xlabel('T AS [%]');
ylabel('Vpp [uV]');
set(gca, 'FontSize', 15);
set(gca, 'XTick', xticks_vec);
set(gca, 'XTickLabel', H_AS/(2*pi*R_nerve*1e+3)*100);
