%% PLOT FIGURES 2, S1
% the data must be simulated with the script run_simulations/sim2_SUAPs_along_fiber.m
% refer to that script for the nomenclature used to store simulated data
DIR_DATA = ""; % TO SET
dir_experiments = fullfile(DIR_DATA, "SUAPs_analysis_along_fiber");

%% PLOT SPATIAL MAPS OF SUAP AMPLITUDES
% NOTE: large-myelinated fiber = 8 um diameter
%       small-myelinated fiber = 4 um diameter
%       unmyelinated fiber     = 1 um diameter

fiber_morpho8 = get_fiber_morpho("ascent_mrg", 8);
fiber_morpho4 = get_fiber_morpho("ascent_mrg", 4);

% uncomment one of the code blocks below to plot the desired figure
% =============================================================================================================================
% INTRAFASCICULAR TIME - LARGE-MYELINATED FIBER
d_fiber = 8;
h_fig = (fiber_morpho8.length_internode/fiber_morpho4.length_internode)*8;
y_el = -[0:2:50];
dz = 5; % um
[z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
levels = [10:10:60];
clim_vec = [11 115];
opts = arrayfun(@(x) {"TIME"        1      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % INTRAFASCICULAR TIME - SMALL-MYELINATED FIBER
% d_fiber = 4;
% h_fig = 8;
% y_el = -[0:2:50];
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [10:10:60];
% clim_vec = [7 67];
% opts = arrayfun(@(x) {"TIME"          1      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % INTRAFASCICULAR TIME - UNMYELINATED FIBER
% d_fiber = 1;
% h_fig = 8;
% y_el = -[0:2:50];
% z_shift = [0; 0.5];
% z_shift_m = z_shift*1e-6;
% levels = [];
% clim_vec = [1.6 28];
% opts = arrayfun(@(x) {"TIME"        1      4       200     0       nan         "ascent_sundt15"         d_fiber         0          z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % INTRAFASCICULAR uNG - LARGE-MYELINATED FIBER
% d_fiber = 8;
% h_fig = (fiber_morpho8.length_internode/fiber_morpho4.length_internode)*8;
% y_el = -[0:2:50];
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [10:10:60];
% clim_vec = [11 115];
% opts = arrayfun(@(x) {"uNG"         2      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % INTRAFASCICULAR uNG - SMALL-MYELINATED FIBER
% d_fiber = 4;
% h_fig = 8;
% y_el = -[0:2:50];
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [10:10:60];
% clim_vec = [7 67];
% opts = arrayfun(@(x) {"uNG"           2      4       200     0       nan         "ascent_mrg"         d_fiber         0          z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % INTRAFASCICULAR uNG - UNMYELINATED FIBER
% d_fiber = 1;
% h_fig = 8;
% y_el = -[0:2:50];
% z_shift = [0; 0.5];
% z_shift_m = z_shift*1e-6;
% levels = [];
% clim_vec = [1.6 28];
% opts = arrayfun(@(x) {"uNG"         2      4       200     0       nan         "ascent_sundt15"         d_fiber         0          z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % EXTRAFASCICULAR TIME - LARGE-MYELINATED FIBER
% d_fiber = 8;
% h_fig = (fiber_morpho8.length_internode/fiber_morpho4.length_internode)*8;
% y_el = -190-[10:2:60];
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [0:2:60];
% clim_vec = [3 15];
% opts = arrayfun(@(x) {"TIME"        1      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % EXTRAFASCICULAR TIME - SMALL-MYELINATED FIBER
% d_fiber = 4;
% h_fig = 8;
% y_el = -190-[10:2:60];
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [3:1:7];
% clim_vec = [2 8];
% opts = arrayfun(@(x) {"TIME"        1      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % EXTRAFASCICULAR TIME - UNMYELINATED FIBER
% d_fiber = 1;
% h_fig = 8;
% y_el = -190-[10:2:60];
% z_shift = [0; 0.5];
% z_shift_m = z_shift*1e-6;
% levels = [];
% clim_vec = [0.4 3];
% opts = arrayfun(@(x) {"TIME"        1      4       200     0       10          "ascent_sundt15"         d_fiber         -190       z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % EXTRAFASCICULAR uNG - LARGE-MYELINATED FIBER
% d_fiber = 8;
% h_fig = (fiber_morpho8.length_internode/fiber_morpho4.length_internode)*8;
% y_el = -190-[10:2:60];
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [0:0.3:2.4];
% clim_vec = [1.4 4];
% opts = arrayfun(@(x) {"uNG"         2      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % EXTRAFASCICULAR uNG - SMALL-MYELINATED FIBER
% d_fiber = 4;
% h_fig = 8;
% y_el = -190-[10:2:60];
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [0.9:0.1:2]; 
% clim_vec = [0.8 1.8];
% opts = arrayfun(@(x) {"uNG"         2      4       200     0       10          "ascent_mrg"         d_fiber         -190       z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % EXTRAFASCICULAR uNG - UNMYELINATED FIBER
% d_fiber = 1;
% h_fig = 8;
% y_el = -190-[10:2:60];
% z_shift = [0; 0.5];
% z_shift_m = z_shift*1e-6;
% levels = [];
% clim_vec = [0.16 0.6];
% opts = arrayfun(@(x) {"uNG"         2      4       200     0       10          "ascent_sundt15"         d_fiber         -190       z_shift     x}, y_el, 'UniformOutput', false);
% =============================================================================================================================
% % EXTRANEURAL CUFF - LARGE-MYELINATED FIBER
% d_fiber = 8;
% h_fig = (fiber_morpho8.length_internode/fiber_morpho4.length_internode)*8;
% d_as_fasc = [0:20:200 300:100:800];
% y_el = -800+d_as_fasc; % NOTE: is the y-position of the fascicle center!
% d_fasc_fib = 200;
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [];
% clim_vec = [5.95 6.5];
% opts = arrayfun(@(x) {"Cortec_Cuff"    2      4       200     x       nan         "ascent_mrg"         d_fiber         x-200+d_fasc_fib       z_shift     nan}, y_el, 'UniformOutput', false);
% y_el = -y_el;
% =============================================================================================================================
% % EXTRANEURAL CUFF - SMALL-MYELINATED FIBER
% d_fiber = 4;
% h_fig = 8;
% d_as_fasc = [0:20:200 300:100:800];
% y_el = -800+d_as_fasc; % NOTE: is the y-position of the fascicle center!
% d_fasc_fib = 200;
% dz = 5; % um
% [z_shift, z_shift_m] = dz_to_zshift("ascent_mrg", d_fiber, dz);
% levels = [];
% clim_vec = [1.86 2.41];
% opts = arrayfun(@(x) {"Cortec_Cuff"    2      4       200     x       nan         "ascent_mrg"         d_fiber         x-200+d_fasc_fib       z_shift     nan}, y_el, 'UniformOutput', false);
% y_el = -y_el;
% =============================================================================================================================
% % EXTRANEURAL CUFF - UNMYELINATED FIBER
% d_fiber = 1;
% h_fig = 8;
% d_as_fasc = [0:20:200 300:100:800];
% y_el = -800+d_as_fasc; % NOTE: is the y-position of the fascicle center!
% d_fasc_fib = 200;
% z_shift = [0; 0.5];
% z_shift_m = z_shift*1e-6;
% levels = [];
% clim_vec = [0.024 0.05];
% opts = arrayfun(@(x) {"Cortec_Cuff"    2      4       200     x       nan         "ascent_sundt15"         d_fiber         x-200+d_fasc_fib       z_shift     nan}, y_el, 'UniformOutput', false);
% y_el = -y_el;
% =============================================================================================================================

N1 = zeros(length(z_shift), length(y_el));
TN1 = zeros(length(z_shift), length(y_el));
P1 = zeros(length(z_shift), length(y_el));
TP1 = zeros(length(z_shift), length(y_el));
Vpp = zeros(length(z_shift), length(y_el));
Tpp = zeros(length(z_shift), length(y_el));
for i=1:length(y_el)
    dir_model = fullfile(dir_experiments, DataHash(opts{i}, 'MD5'));

    if strcmp(opts{i}{1},'Cortec_Cuff')
        SUAP_t = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP_t;
        SUAP1 = load(fullfile(dir_model,"SUAP_1_1.mat")).SUAP;
        SUAP2 = load(fullfile(dir_model,"SUAP_1_2.mat")).SUAP;
        SUAP3 = load(fullfile(dir_model,"SUAP_1_3.mat")).SUAP;
    else
        load(fullfile(dir_model,"SUAP_1_1.mat"));
    end

    for j=1:length(z_shift)
        if strcmp(opts{i}{1},'Cortec_Cuff')
            SUAP_j = (2*SUAP2{j}-SUAP1{j}-SUAP3{j})/2;
        else
            SUAP_j = SUAP{j};
        end
        
        % Vpp(j, i) = max(SUAP_j)-min(SUAP_j);

        [N1(j,i), TN1_idx] = min(SUAP_j);
        TN1(j,i) = SUAP_t{j}(TN1_idx);
        [P1(j,i), TP1_idx] = max(SUAP_j(SUAP_t{j}<TN1(j,i)));
        TP1(j,i) = SUAP_t{j}(TP1_idx);
        Vpp(j,i) = P1(j,i)-N1(j,i);
        Tpp(j,i) = TN1(j,i)-TP1(j,i);
    end
end

figure;
[X,Y] = meshgrid(-y_el, z_shift_m*1e+6);
surf(X, Y, Vpp, 'LineStyle', 'none', 'FaceColor', 'interp');
hold on;
if ~isempty(levels)
    [C,h] = contour(X, Y, Vpp, levels, 'w'); % 'ShowText','on'
    clabel(C,h,'FontSize',15,'Color','w')
    h_obj = findall(gcf, 'type', 'Contour');
    h_obj.ZLocation = 1000;
    h_obj.LineWidth = 1;
end
h_cb = colorbar;
h_cb.Title.String = 'Vpp [uV]';
if ~isempty(clim_vec)
    clim(clim_vec);
end
colormap jet
xlabel('Distance from fiber [um]', 'FontSize', 15);
ylabel('Internode position [um]', 'FontSize', 15);
set(gca, 'color', 'none', 'FontSize', 15);
ylim(z_shift_m([1 end])*1e+6);
view(2);
grid off;
box on;
set(gcf, 'units', 'centimeters');
pos = get(gcf,'Position');
set(gcf, 'Position', [pos(1) pos(2)-(h_fig-pos(4)) pos(3) h_fig]);
