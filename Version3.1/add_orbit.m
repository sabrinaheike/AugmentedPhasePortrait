function[ ] = add_orbit(fin2,gin2,Xo,Yo,iterations,options)
    arguments 
        fin2 (1,1) function_handle
        gin2 (1,1) function_handle
        Xo (1,1) {mustBeNumeric, mustBeReal}
        Yo (1,1) {mustBeNumeric, mustBeReal}
        iterations (1,1) {mustBeNumeric, mustBeReal}
        options.OrbitAxes (1,1)
    end

    if isfield(options,'OrbitAxes')
        ax = options.OrbitAxes;
    end

    if isfield(options,'OrbitAxes')
        hold(ax,"on")
        plot(ax,Xo,Yo,'-*','Color','m',"MarkerSize",20,"HandleVisibility","off")
    else
        hold on
        plot(Xo,Yo,'-*','Color','m',"MarkerSize",20,"HandleVisibility","off")
    end
    a = ones(iterations+1,2);
    a(1,1) = Xo;
    a(1,2) = Yo;
    for c = 2:iterations+1
        a(c,1) = fin2(a(c-1,1),(a(c-1,2)));
        a(c,2) = gin2(a(c-1,1),(a(c-1,2)));

    end
    if isfield(options,'OrbitAxes')
        plot(ax, a(:,1),a(:,2),'-','Color','m','Marker','o',"LineWidth",1,"DisplayName","orbit","HandleVisibility","off")
        hold(ax,"off")
    else
        plot(a(:,1),a(:,2),'-','Color','m','Marker','o',"LineWidth",1,"DisplayName","orbit","HandleVisibility","off")
        hold(ax,"off")
    end
end
