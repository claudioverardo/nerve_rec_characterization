function colors = data2color(CData, CLim, colormap)
% DATA2COLOR Map a set of values (CData) in a set of colors according to two
% correspondent ranges of data (CLim) and colors (colormap).
% 
%   color = data2color(CData, CLim, colormap)
% 
% https://it.mathworks.com/help/matlab/creating_plots/change-mapping-of-data-values-into-the-colormap.html

    n_colors = size(colormap,1);
    color_idcs = round((CData-CLim(1))/diff(CLim)*(n_colors-1))+1;
    
    color_idcs(color_idcs<1) = 1;
    color_idcs(color_idcs>n_colors) = n_colors;
    
    colors = colormap(color_idcs,:);

end