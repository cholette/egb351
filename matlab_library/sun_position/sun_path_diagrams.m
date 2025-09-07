function [fig,ax] = sun_path_diagrams(T)

    az = T.SolarAzimuthAngle;
    z = T.SolarZenithAngle;

    fig = figure(name='Sun paths');
    tiledlayout(2,1)
    ax1 = nexttile;
    
    doys = 0:5:365;
    C = hsv(length(doys)); % circular colormap
    for ii = 1:length(doys)
        d = doys(ii);
        plot(az(T.doy==d),90-z(T.doy==d),'-.','Color',C(ii,:));
        hold on
    end
    zlabel('Solar azimuth')
    ylabel('Solar elevation angle (deg)')
    title('(+ is West of South)')
    ylim([0,90])
    colormap('hsv');
    c = colorbar;
    clim([min(doys) max(doys)])
    c.Label.String = 'Day of year';
    hold off
    
    ax2 = nexttile;
    azn = az_zero_north(az);
    for ii = 1:length(doys)
        d = doys(ii);
        plot(azn(T.doy==d),90-z(T.doy==d),'-.','Color',C(ii,:));
        hold on
    end
    ylabel('Solar elevation angle (deg)')
    xlabel('Solar azimuth')
    title('(+ is West of North)')
    ylim([0,90])
    colormap('hsv');
    c = colorbar;
    clim([min(doys) max(doys)])
    c.Label.String = 'Day of year';
    hold off

    ax = [ax1,ax2];
    linkaxes(ax,'x')