function [] = plot_heatmap(map,range,colmap,maptitle,cluster_count, to_save)
%plot and save heatmaps
%   plot values within the range with the specified color map 
%   with title
%   denote cluster boundary using cluster count
%   save to 'figures/' if to_save=true(1)
    imagesc(map);
	colormap(colmap);
	caxis(range);
	colorbar;
	title(maptitle, 'interpreter', 'none');
    axis square
	xticks(cumsum(cluster_count));
	yticks(cumsum(cluster_count));
	grid on
    drawnow;
    filename=strrep(maptitle,{' '}, '_');
    fn = char(strcat('figures/', filename{1}));
    if to_save
        saveas(gcf, fn, 'jpg');
    end

end

