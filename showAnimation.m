function [] = showAnimation(coords_t, L)
        nPart = size(coords_t,2);
        c = linspace(0,255,nPart);
        scatter3 ( coords_t(1,:), coords_t(2,:), coords_t(3,:), [], c, 'filled')
        %title ( sprintf ( 'Step %d\n', step ) );
        xlim([0 L]);
        ylim([0 L]);
        zlim([0 L]);
        view(30,40);
        pause (0.05);
end