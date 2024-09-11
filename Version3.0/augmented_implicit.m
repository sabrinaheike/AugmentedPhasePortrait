% Implicit Augmented phase portrait
% Input: (f,g,minx, maxx, miny, maxy)

% f = @(x,y) f(x,y); right-hand side of the X-equation f = f(x,y)
% g = @(x,y) g(x,y); right-hand side of the Y-equation g = g(x,y) 
% minx = number; minimum x-value to be considered for plotting
% maxx = number; maximum x-value to be considered for plotting
% miny = number; minimum y-value to be considered for plotting
% maxy = number; maximum y-value to be considered for plotting

% optional parameters: add as Name-Value pairs
% Accuracy = number; optional parameter for number of arrows and signs of 
%   next-iterate operators; default value = 20
% Cutoff Values = 1x2 double; if numbers given then dash-dotted curves are
%   added in plot to show points (x,y) such that f(x,y)=cutoffx and
%   g(x,y)=cutoffy
% PortraitAxes = axes object; specifies the axes to plot the augmented
%   phase portrait
% NullclineRange = 1x4 double; Optional range over which to calculate the
%   nullclines. The root curves and symbols of next iterate operators
%   correspond only to the sections of nullclines that are plotted. In order
%   to make the root curves and symbols extend further, add a nullcline calculation 
%   range which is larger than the minimum and maximum values for plotting.
%   The nullcline calculation values are to be inputed in the format
%   [minx, maxx, miny, maxy].

function [] = augmented_implicit(fin,gin, minx, maxx, miny, maxy, options)
    
    % Argument validation
    arguments
        fin (1,1) function_handle
        gin (1,1) function_handle
        minx (1,1) {mustBeNumeric, mustBeReal}
        maxx (1,1) {mustBeNumeric, mustBeReal}
        miny (1,1) {mustBeNumeric, mustBeReal}
        maxy (1,1) {mustBeNumeric, mustBeReal}
        options.Accuracy (1,1) {mustBeNumeric, mustBeReal} = 20
        options.CutoffValues (1,2) double {mustBeNumeric,mustBeReal}
        options.PortraitAxes (1,1) = axes()
        options.NullclineRange (1,4) double {mustBeNumeric, mustBeReal}
    end

    % Check if optional parameters are provided

    % Nullcline calculation range
        % the variables with an underscore (ex. max_y) are for calculating
        % the nulllclines. The values without (ex. maxy) are used for
        % plotting.
    if isfield(options,'NullclineRange') % if a nullcline range is given
        min_x = options.NullclineRange(1); % use values given by user
        max_x = options.NullclineRange(2);
        min_y = options.NullclineRange(3);
        max_y = options.NullclineRange(4);
    else % use the same values as the plotting range
        min_x = minx;
        max_x = maxx;
        min_y = miny;
        max_y = maxy;
    end
    
    % Set accuracy and axes to the argument given, if no values are given
    % then the default values are used
    acc = options.Accuracy;
    ax = options.PortraitAxes;

    % Create nullcline equations
    Xfunction = @(x,y) x;
    Yfunction = @(x,y) y;
    
    Xnullclines = @(x,y) fin(x,y) - Xfunction(x,y);
    Ynullclines = @(x,y) gin(x,y) - Yfunction(x,y);

    % Find nullclines using fimplicit
    PlotTolerance = 0.01; % fimplicit looks for nullclines in a slightly larger range than the one inputted so that nullclines on the edges of the range are not missed. 
    XnullclinesPlot = fimplicit(ax,Xnullclines,[min_x-PlotTolerance max_x+PlotTolerance min_y-PlotTolerance max_y+PlotTolerance], 'LineStyle', 'none', 'HandleVisibility','off','MeshDensity',300);
    hold(ax,"on")
    YnullclinesPlot = fimplicit(ax,Ynullclines,[min_x-PlotTolerance max_x+PlotTolerance min_y-PlotTolerance max_y+PlotTolerance], 'LineStyle', 'none', 'HandleVisibility','off','MeshDensity',300);

    % Extract X and Y data
    XnullclinesXdata = round(XnullclinesPlot.XData,5);
    XnullclinesYdata = round(XnullclinesPlot.YData,5);
    YnullclinesXdata = round(YnullclinesPlot.XData,5);
    YnullclinesYdata = round(YnullclinesPlot.YData,5);

    % Separate Nullcline Data Into Branches
        % a) X nullclines
    
    diffx = diff([XnullclinesXdata; XnullclinesYdata],1,2); % diff between x vals in row 1 and diff between y vals in row 2
    diffx(diffx>0) = 1; % set positive differences to 1
    diffx(diffx<0) = -1; % set negative differences to -1
    recordBreaksX = zeros(length(XnullclinesXdata),1); % Value of 1 signifies the start of a new nullcline, value of 2 signifies a point to be removed 
    variable = []; % if variable = 1 the nullcline is labelled as a function of x, 0 for a function of y 
    for i = 1:size(diffx,2)
        if isnan(diffx(1,i)) || isnan(diffx(2,i)) % NaN value
            recordBreaksX(i+1) = 2; % remove the NaN value
            recordBreaksX(i+2) = 1; % point after NaN starts a new nullcline branch
            i = i + 1; % no need to check next difference, since it will also be NaN
            variable = [];
        elseif diffx(1,i) == 0 && diffx(2,i) == 0 %same value
            recordBreaksX(i+1) = 2; % remove value
        elseif diffx(1,i) == 0 % constant x values
            if variable == 1 % if nullcline is a function of x
                recordBreaksX(i+1) = 1; % mark the next value as the start of a new nulcline
                variable = []; % reset the variable
            else % mark the nullcline as a function of y
                variable = 0;
            end
        elseif diffx(2,i) == 0 % constant y values
            if variable == 0 % if nullcline is a function of y
                recordBreaksX(i+1) = 1; % mark the next value as the start of a new nulcline
                variable = []; % reset the variable
            else % mark the nullcline as a function of x 
                variable = 1;
            end
        elseif i>=2 && recordBreaksX(i)~=1 % not the start of a new nullcline
            firstdiff_x = diffx(1,i);
            firstdiff_y = diffx(2,i);
            if recordBreaksX(i) == 2
                if recordBreaksX(i-1)==1
                    seconddiff_x = firstdiff_x;
                    seconddiff_y = firstdiff_y;
                else
                    seconddiff_x = diffx(1,i-2);
                    seconddiff_y = diffx(2,i-2);
                end
            else
                seconddiff_x = diffx(1,i-1);
                seconddiff_y = diffx(2,i-1);
            end
            if firstdiff_x - seconddiff_x ~= 0 && firstdiff_y - seconddiff_y ~= 0 % both x and y change direction
                recordBreaksX(i+1) = 1;
                variable = [];
            elseif firstdiff_x - seconddiff_x ~= 0 % x values change direction
                if variable == 1
                    recordBreaksX(i+1) = 1;
                    variable = [];
                else
                    variable = 0;
                end
            elseif firstdiff_y - seconddiff_y ~= 0 % y values change direction
                if variable == 0
                    recordBreaksX(i+1) = 1;
                    variable = [];
                else
                    variable = 1;
                end
            end
        end
    end

    % remove duplicate values
    XnullclinesXdata = XnullclinesXdata(recordBreaksX~=2);
    XnullclinesYdata = XnullclinesYdata(recordBreaksX~=2);

    recordBreaksX = recordBreaksX(recordBreaksX~=2);

    % number of x nullclines
    if ~isempty(XnullclinesXdata)
        nr_nullclinesX = sum(recordBreaksX == 1) + 1;
    else 
        nr_nullclinesX = 0;
    end

    % record indices of nullcline breaks 
    Xidx = find(recordBreaksX==1);
    if ~isempty(XnullclinesXdata)
        Xidx = [1; Xidx; length(XnullclinesXdata)+1];
    end

        % create array to store data in separate branches
    XbranchesXdata = nan(2*length(XnullclinesXdata),nr_nullclinesX);
    XbranchesYdata = nan(2*length(XnullclinesXdata),nr_nullclinesX);

        % sort the data into the array 
    for i = 2:length(Xidx)
        nullclineValuesX = XnullclinesXdata(Xidx(i-1):Xidx(i)-1);
        nullclineValuesY = XnullclinesYdata(Xidx(i-1):Xidx(i)-1);
        XbranchesXdata(1:length(nullclineValuesX),i-1) = nullclineValuesX;
        XbranchesYdata(1:length(nullclineValuesX),i-1) = nullclineValuesY;
    end

        % b) Y nullclines

    diffy = diff([YnullclinesXdata; YnullclinesYdata],1,2);
    diffy(diffy>0) = 1;
    diffy(diffy<0) = -1;
    recordBreaksY = zeros(length(YnullclinesXdata),1);
    variable = [];
    for i = 1:size(diffy,2)
        if isnan(diffy(1,i)) || isnan(diffy(2,i))
            recordBreaksY(i+1) = 2;
            recordBreaksY(i+2) = 1;
            variable = [];
        elseif diffy(1,i) == 0 && diffy(2,i) == 0 %same value
            recordBreaksY(i+1) = 2; % remove value
        elseif diffy(1,i) == 0
            if variable == 1
                recordBreaksY(i+1) = 1;
                variable = [];
            else
                variable = 0;
            end
        elseif diffy(2,i) == 0
            if variable == 0
                recordBreaksY(i+1) = 1;
                variable = [];
            else
                variable = 1;
            end
        elseif i>=2 && recordBreaksY(i)~=1
            firstdiff_x = diffy(1,i);
            firstdiff_y = diffy(2,i);
            if recordBreaksY(i) == 2
                if recordBreaksY(i-1)==1
                    seconddiff_x = firstdiff_x;
                    seconddiff_y = firstdiff_y;
                else
                    seconddiff_x = diffy(1,i-2);
                    seconddiff_y = diffy(2,i-2);
                end
            else
                seconddiff_x = diffy(1,i-1);
                seconddiff_y = diffy(2,i-1);
            end
            if firstdiff_x - seconddiff_x ~= 0 && firstdiff_y - seconddiff_y ~= 0 % both x and y change direction
                recordBreaksY(i+1) = 1;
                variable = [];
            elseif firstdiff_x - seconddiff_x ~= 0
                if variable == 1
                    recordBreaksY(i+1) = 1;
                    variable = [];
                else
                    variable = 0;
                end
            elseif firstdiff_y - seconddiff_y ~= 0
                if variable == 0
                    recordBreaksY(i+1) = 1;
                    variable = [];
                else
                    variable = 1;
                end
            end
        end
    end

    % remove duplicate values
    YnullclinesXdata = YnullclinesXdata(recordBreaksY~=2);
    YnullclinesYdata = YnullclinesYdata(recordBreaksY~=2);

    recordBreaksY = recordBreaksY(recordBreaksY~=2);

    % number of x nullclines
    if ~isempty(YnullclinesXdata)
        nr_nullclinesY = sum(recordBreaksY == 1) + 1;
    else 
        nr_nullclinesY = 0;
    end

    % record indices of nullcline breaks 
    Yidx = find(recordBreaksY==1);
    if ~isempty(YnullclinesXdata)
        Yidx = [1; Yidx; length(YnullclinesXdata)+1];
    end

    YbranchesXdata = nan(length(YnullclinesXdata),nr_nullclinesY);
    YbranchesYdata = nan(length(YnullclinesXdata),nr_nullclinesY);

    for i = 2:length(Yidx)
        nullclineValuesX = YnullclinesXdata(Yidx(i-1):Yidx(i)-1);
        nullclineValuesY = YnullclinesYdata(Yidx(i-1):Yidx(i)-1);
        YbranchesXdata(1:length(nullclineValuesX),i-1) = nullclineValuesX;
        YbranchesYdata(1:length(nullclineValuesX),i-1) = nullclineValuesY;
    end


    % Plot nullcline branches 
        % a) X nullclines
        % Set colours
    if nr_nullclinesX <=3
        XColours = [0.35, 0.5, 1; 0, 0, 1; 0, 0, 0.6];
    else 
        XColours = colormap(winter(size(XnullclineInfo,1)+4));
    end
    
        % Plot
    for i = 1:nr_nullclinesX
        if ~isnan(XbranchesXdata(1,i))
            colr = XColours(i,:);
            plot(ax,XbranchesXdata(:,i),XbranchesYdata(:,i),'--', 'Color', colr, 'LineWidth', 3, 'DisplayName',strcat('X-Nullcline',num2str(i)))
        end
    end
        % b) Y nullclines
        % Set colours
    if nr_nullclinesY <=3
        YColours = [0.75, 0, 0; 1, 0, 0; 0.5, 0, 0];
    else
        YColours = colormap(autumn(size(YnullclineInfo,1)+4));
    end

    for i = 1:nr_nullclinesY
        if ~isnan(YbranchesXdata(1,i))
            colr = YColours(i,:);
            plot(ax,YbranchesXdata(:,i),YbranchesYdata(:,i),'--', 'Color', colr, 'LineWidth', 3, 'DisplayName',strcat('Y-Nullcline',num2str(i)))
        end
    end

    % Construct Grid
    xaxis = minx:(maxx-minx)/acc:maxx; 
    yaxis = miny:(maxy-miny)/acc:maxy;
    nr = max(nr_nullclinesX,nr_nullclinesY);
   
    % Draw Direction Field
    xaxish = xaxis(1:nr+2:end);
    yaxish = yaxis(1:nr+2:end);

    drawArrow = @(w,z, varargin) quiver(ax, w(1),z(1),w(2)-w(1),z(2)-z(1),2, 'MaxHeadSize', 12, varargin{:});    
    
    for md=1:length(xaxish)
        for nd = 1:length(yaxish)
            % points to be evaluated
            xt = xaxish(md);
            yt = yaxish(nd);
            fval = fin(xt,yt);
            gval = gin(xt,yt);
            % end of the arrowlength to be evaluated (to avoid arrows
            % pointing across nullclines in the wrong direction)
            arrowlex = (maxx-minx)/(2*acc); 
            arrowley = (maxy-miny)/(2*acc);
            xt2 = xt+sign(fval-xt)*2*arrowlex;
            yt2 = yt+sign(gval-yt)*2*arrowley; 
            fval2 = fin(xt2,yt);
            gval2 = gin(xt,yt2);
            % arrowlength
            if xt2 >= minx && xt2 <= maxx && yt2 >= miny && yt2 <= maxy
                if fval>xt && fval2>xt2 && gval>yt  && gval2>yt2
                    hold(ax,"on")
                    drawArrow([xt, xt+arrowlex],[yt, yt],'linewidth',1,'color','k','HandleVisibility','off');     
                    hold(ax,"on")
                    drawArrow([xt, xt],[yt, yt+arrowley], 'linewidth',1,'color','k','HandleVisibility','off');
                elseif fval>xt && fval2>xt2 && gval<yt  && gval2<yt2
                    hold(ax,"on")
                    drawArrow([xt, xt],[yt, yt-arrowley],'linewidth',1,'color','k','HandleVisibility','off');     
                    hold(ax,"on")
                    drawArrow([xt, xt+arrowlex],[yt, yt], 'linewidth',1,'color','k','HandleVisibility','off');
    
                elseif fval<xt  && fval2< xt2 && gval<yt  && gval2<yt2
                    hold(ax,"on")
                    drawArrow([xt, xt],[yt, yt-arrowley],'linewidth',1,'color','k','HandleVisibility','off');     
                    hold(ax,"on")
                    drawArrow([xt, xt-arrowlex],[yt, yt],'linewidth',1,'color','k','HandleVisibility','off');     
                elseif fval<xt  && fval2<xt2 && gval>yt  && gval2>yt2
                    hold(ax,"on")
                    drawArrow([xt, xt],[yt, yt+arrowley],'linewidth',1,'color','k','HandleVisibility','off');     
                    hold(ax,"on")
                    drawArrow([xt, xt-arrowlex],[yt, yt],'linewidth',1,'color','k','HandleVisibility','off');     
                end
            end
        end
    end


    % Construct next iterate operator
        % a) X nullclines
    for i = 1:nr_nullclinesX
        colr = XColours(i,:); %set colour
        
        % x and y data from nullcline
        xvals = XbranchesXdata(:,i);
        xreal = xvals(~isnan(xvals));
        yvals = XbranchesYdata(:,i);
        yreal = yvals(~isnan(yvals));

        % Check if branch is a function of x or y 
        xdir = diff(xreal);
        ydir = diff(yreal);
        xfunc = 0;
        yfunc = 0;
        if all(xdir < 0) || all(xdir>0) % x values change in constant direction
            xfunc = 1; % label as x function
        elseif all(ydir < 0) || all(ydir>0) % y values change in constant direction
            yfunc = 1; % label as y function
        end
        
        % create NIO function
        if xfunc == 1
            nullclineFunction = @(xTest) interp1(xreal,yreal,xTest);
            subNullclineF = @(x,y) nullclineFunction(fin(x,y));
            NIO = @(x,y) gin(x,y) - subNullclineF(x,y);
        elseif yfunc == 1
            nullclineFunction = @(yTest) interp1(yreal,xreal,yTest);
            subNullclineG = @(x,y) nullclineFunction(gin(x,y));
            NIO = @(x,y) fin(x,y) - subNullclineG(x,y);
        end
        
        % plot root curve
        fimplicit(ax,NIO,[minx-PlotTolerance maxx+PlotTolerance miny-PlotTolerance maxy+PlotTolerance],'Color',colr,'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
        
        % evaluate NIO on grid
        if xfunc == 1
            for j = i+1:nr+2:length(xaxis)
                for k = i+2:nr+2:length(yaxis)
                    symbolx = xaxis(j);
                    symboly = yaxis(k);
                    NIOvalue = NIO(symbolx,symboly);
                    if isreal(NIOvalue) && NIOvalue>0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','^','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    elseif isreal(NIOvalue) && NIOvalue<0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','v','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    else
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',10,'Marker','pentagram','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    end
                end
            end
        elseif yfunc ==1
           for j = i+1:nr+2:length(xaxis)
                for k = i+2:nr+2:length(yaxis)
                    symbolx = xaxis(j);
                    symboly = yaxis(k);
                    NIOvalue = NIO(symbolx,symboly);
                    if isreal(NIOvalue) && NIOvalue>0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','>','LineStyle','-','LineWidth',0.5,'HandleVisibility','off');
                    elseif isreal(NIOvalue) && NIOvalue<0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','<','LineStyle','-','LineWidth',0.5,'HandleVisibility','off');
                    else
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',10,'Marker','pentagram','LineStyle','-','LineWidth',0.5,'HandleVisibility','off');
                    end
                end
            end 
        end
    end

        % b) Y nullclines
    for i = 1:nr_nullclinesY
        colr = YColours(i,:);
        
        % x and y data from nullcline
        xvals = YbranchesXdata(:,i);
        xreal = xvals(~isnan(xvals));
        yvals = YbranchesYdata(:,i);
        yreal = yvals(~isnan(yvals));

        % Check if branch is a function of x or y 
        xdir = diff(xreal);
        ydir = diff(yreal);
        xfunc = 0;
        yfunc = 0;
        if all(xdir < 0) || all(xdir>0)
            xfunc = 1;
        elseif all(ydir < 0) || all(ydir>0)
            yfunc = 1;
        end
        
        % create NIO function
        if xfunc == 1
            nullclineFunction = @(xTest) interp1(xreal,yreal,xTest);
            subNullclineF = @(x,y) nullclineFunction(fin(x,y));
            NIO = @(x,y) gin(x,y) - subNullclineF(x,y);
        elseif yfunc == 1
            nullclineFunction = @(yTest) interp1(yreal,xreal,yTest);
            subNullclineG = @(x,y) nullclineFunction(gin(x,y));
            NIO = @(x,y) fin(x,y) - subNullclineG(x,y);
        end

        % plot root curve
        fimplicit(ax,NIO,[minx-PlotTolerance maxx+PlotTolerance miny-PlotTolerance maxy+PlotTolerance],'Color',colr,'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
        
        % evaluate NIO on grid
        if xfunc ==1
            for j = i+2:nr+2:length(xaxis)
                for k = i+1:nr+2:length(yaxis)
                    symbolx = xaxis(j);
                    symboly = yaxis(k);
                    NIOvalue = NIO(symbolx,symboly);
                    if isreal(NIOvalue) && NIOvalue>0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','^','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    elseif isreal(NIOvalue) && NIOvalue<0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','v','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    else
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',10,'Marker','pentagram','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    end
                end
            end
        elseif yfunc == 1
            for j = i+2:nr+2:length(xaxis)
                for k = i+1:nr+2:length(yaxis)
                    symbolx = xaxis(j);
                    symboly = yaxis(k);
                    NIOvalue = NIO(symbolx,symboly);
                    if isreal(NIOvalue) && NIOvalue>0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','>','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    elseif isreal(NIOvalue) && NIOvalue<0
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','<','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    else
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',10,'Marker','pentagram','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    end
                end
            end
        end
    end

    % Cutoff functions
    if isfield(options,'CutoffValues')
        cutoffvalx = @(x,y) options.CutoffValues(1); % function handle containing cutoff value
        fcutf = @(x,y) fin(x,y) - cutoffvalx(x,y);
        fimplicit(ax,fcutf,':', 'Color',[0.45 0.85 0.45], 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('f=',num2str(options.CutoffValues(1))))
        cutoffvaly = @(x,y) options.CutoffValues(2);
        gcutf = @(x,y) gin(x,y) - cutoffvaly(x,y);
        fimplicit(ax,gcutf,':', 'Color',[0 0.45 0], 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('g=',num2str(options.CutoffValues(2))))

    end

    % Lines along coordinate axes
    xline(ax,0,"LineStyle",'--',"Color",[0.5 0.5 0.5],'HandleVisibility','off')
    yline(ax,0,"LineStyle",'--',"Color",[0.5 0.5 0.5],'HandleVisibility','off')
    
    % Set axes limits
    ax.XLim = [minx maxx];
    ax.YLim = [miny maxy];

    % Title, axes labels and legend
    xlabel(ax,'X');
    ylabel(ax,'Y', 'rotation', 0);
    title(ax,'Phase Portrait')

    legend(ax)
    set(ax,'FontSize',18)



    