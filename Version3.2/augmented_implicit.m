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
% CombineNullcline = "on" or "off". Optional input to control the automatic
%   combining of nullcline segments which can together form 1 continuous
%   nullcline function. The default setting is "on"

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
        options.CombineNullcline (1,1) string = "on"
        options.OptimizeNullclineCombination (1,1) string = "off"
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
    mesh = 100;
    PlotToleranceX = (max_x-min_x)/1000; % fimplicit looks for nullclines in a slightly larger range than the one inputted so that nullclines on the edges of the range are not missed. 
    PlotToleranceY = (max_y-min_y)/1000;
    XnullclinesPlot = fimplicit(ax,Xnullclines,[min_x-PlotToleranceX max_x+PlotToleranceX min_y-PlotToleranceY max_y+PlotToleranceY], 'LineStyle', 'none', 'HandleVisibility','off','MeshDensity',mesh);
    hold(ax,"on")
    YnullclinesPlot = fimplicit(ax,Ynullclines,[min_x-PlotToleranceX max_x+PlotToleranceX min_y-PlotToleranceY max_y+PlotToleranceY], 'LineStyle', 'none', 'HandleVisibility','off','MeshDensity',mesh);

    % Extract X and Y data
    
    XnullclinesXdata = XnullclinesPlot.XData;
    XnullclinesYdata = XnullclinesPlot.YData;
    YnullclinesXdata = YnullclinesPlot.XData;
    YnullclinesYdata = YnullclinesPlot.YData;

    % Round points for seperation
    range = (max_x - min_x) + (max_y - min_y);
    rangeScale = round(log10(range));
    decimal = -rangeScale + 4;

    XnullclinesXdataRounded = round(XnullclinesXdata,decimal);
    XnullclinesYdataRounded = round(XnullclinesYdata,decimal);
    YnullclinesXdataRounded = round(YnullclinesXdata,decimal);
    YnullclinesYdataRounded = round(YnullclinesYdata,decimal);

    % Separate Nullcline Data Into Branches
    recordBreaksX = zeros(length(XnullclinesXdata),1); % Value of 1 signifies the start of a new nullcline, value of 2 signifies a point to be removed         

    recordBreaksX = splitNullclineData(XnullclinesXdata, XnullclinesYdata, recordBreaksX);

    recordBreaksX = splitNullclineData(XnullclinesXdataRounded, XnullclinesYdataRounded, recordBreaksX);
    % remove duplicate values
    XnullclinesXdata = XnullclinesXdata(recordBreaksX~=3);
    XnullclinesYdata = XnullclinesYdata(recordBreaksX~=3);

    XnullclinesXdataRounded = XnullclinesXdataRounded(recordBreaksX~=3);
    XnullclinesYdataRounded = XnullclinesYdataRounded(recordBreaksX~=3);

    recordBreaksX = recordBreaksX(recordBreaksX~=3);

    % number of x nullclines
    if ~isempty(XnullclinesXdata)
        % number of nullclines is the number of breaks due to NaN values plus the number of breaks due to changing direction plus 1
        nr_nullclinesX = sum(recordBreaksX == 1) + sum(recordBreaksX == 2) + 1;
    else 
        nr_nullclinesX = 0;
    end

    % record indices of nullcline breaks 
    Xidx = [find(recordBreaksX==1); find(recordBreaksX==2)];
    Xidx = sort(Xidx);
    if ~isempty(XnullclinesXdata)
        Xidx = [1; Xidx; length(XnullclinesXdata)+1]; % indices indicate when new nullcline branches start, plus 1 additional index at the end (would be the start of the next nullcline if there was one)
    end

        % create array to store data in separate branches
    XbranchesXdata = nan(nr_nullclinesX*length(XnullclinesXdata),nr_nullclinesX);
    XbranchesYdata = nan(nr_nullclinesX*length(XnullclinesXdata),nr_nullclinesX);

    XbranchesXdataRounded = nan(nr_nullclinesX*length(XnullclinesXdata),nr_nullclinesX);
    XbranchesYdataRounded = nan(nr_nullclinesX*length(XnullclinesXdata),nr_nullclinesX);

        % sort the data into the array 
    shortestnullcline = 10; % don't want to do calculations for nullcline segments which are too short to properly see
    varianceX = (maxx-minx)/100; % how far we allow the x-values to change in the next iteration from a nullcline point
    for i = 2:length(Xidx)
        if recordBreaksX(Xidx(i-1)) == 2 || Xidx(i-1) == 1 % don't share a common point between nullcline branches
            nullclineValuesX = XnullclinesXdata(Xidx(i-1):Xidx(i)-1);
            nullclineValuesY = XnullclinesYdata(Xidx(i-1):Xidx(i)-1);

            nullclineValuesXrounded = XnullclinesXdataRounded(Xidx(i-1):Xidx(i)-1);
            nullclineValuesYrounded = XnullclinesYdataRounded(Xidx(i-1):Xidx(i)-1);

        elseif recordBreaksX(Xidx(i-1)) == 1 % can be made to share a common point between nullcline branches, new branch starts due to change in direction not becuase of a NaN
            nullclineValuesX = XnullclinesXdata(Xidx(i-1):Xidx(i)-1);
            nullclineValuesY = XnullclinesYdata(Xidx(i-1):Xidx(i)-1);

            nullclineValuesXrounded = XnullclinesXdataRounded(Xidx(i-1):Xidx(i)-1);
            nullclineValuesYrounded = XnullclinesYdataRounded(Xidx(i-1):Xidx(i)-1);
        end
        if length(nullclineValuesX) > shortestnullcline % if the nullcline is long enough
            
            % for each point on the nullcline, take the absolute difference
            % between the x value of the point, and the x value of the
            % point that it is mapped to. If it is a true nullcline this
            % should be very close to 0 (within our threshold of variance)
            nullValues = fin(nullclineValuesX,nullclineValuesY) - nullclineValuesX; 
            nr_nullValues = abs(nullValues) < varianceX; % indices of true nullcline values
            
            nullclineValuesX = nullclineValuesX(nr_nullValues);
            nullclineValuesY = nullclineValuesY(nr_nullValues);
            
            nullclineValuesXrounded = nullclineValuesXrounded(nr_nullValues);
            nullclineValuesYrounded = nullclineValuesYrounded(nr_nullValues);

            if length(nullclineValuesX) > shortestnullcline %if the nullcline with erroneous points removed is still long enough
                % add values to array 
                XbranchesXdata(1:length(nullclineValuesX),i-1) = nullclineValuesX;
                XbranchesYdata(1:length(nullclineValuesX),i-1) = nullclineValuesY;

                XbranchesXdataRounded(1:length(nullclineValuesX),i-1) = nullclineValuesXrounded;
                XbranchesYdataRounded(1:length(nullclineValuesX),i-1) = nullclineValuesYrounded;
            end
        end
    end

    % update the number of X-nullclines since it might have changed when we
    % removed nullclines that were too short or were not true nullclines
    if ~isempty(XbranchesXdata)
        idx = ~(isnan(XbranchesXdata(1,:)));
        XbranchesXdata = XbranchesXdata(:,idx);
        XbranchesYdata = XbranchesYdata(:,idx);

        XbranchesXdataRounded = XbranchesXdataRounded(:,idx);
        XbranchesYdataRounded = XbranchesYdataRounded(:,idx);

        nr_nullclinesX = sum(idx);
    end

        % b) Y nullclines
    
    recordBreaksY = zeros(length(YnullclinesXdata),1);
    recordBreaksY = splitNullclineData(YnullclinesXdata, YnullclinesYdata, recordBreaksY);

    recordBreaksY = splitNullclineData(YnullclinesXdataRounded, YnullclinesYdataRounded, recordBreaksY);
    
    % remove duplicate values
    YnullclinesXdata = YnullclinesXdata(recordBreaksY~=3);
    YnullclinesYdata = YnullclinesYdata(recordBreaksY~=3);

    YnullclinesXdataRounded = YnullclinesXdataRounded(recordBreaksY~=3);
    YnullclinesYdataRounded = YnullclinesYdataRounded(recordBreaksY~=3);

    recordBreaksY = recordBreaksY(recordBreaksY~=3);

    % number of x nullclines
    if ~isempty(YnullclinesXdata)
        nr_nullclinesY = sum(recordBreaksY == 1) + sum(recordBreaksY == 2) + 1;
    else 
        nr_nullclinesY = 0;
    end

    % record indices of nullcline breaks 
    Yidx = [find(recordBreaksY==1); find(recordBreaksY==2)];
    Yidx = sort(Yidx);
    if ~isempty(YnullclinesXdata)
        Yidx = [1; Yidx; length(YnullclinesXdata)+1];
    end

    YbranchesXdata = nan(nr_nullclinesY*length(YnullclinesXdata),nr_nullclinesY);
    YbranchesYdata = nan(nr_nullclinesY*length(YnullclinesXdata),nr_nullclinesY);

    YbranchesXdataRounded = nan(nr_nullclinesY*length(YnullclinesXdata),nr_nullclinesY);
    YbranchesYdataRounded = nan(nr_nullclinesY*length(YnullclinesXdata),nr_nullclinesY);

    % sort the data into the array
    varianceY = (maxy-miny)/100;
    for i = 2:length(Yidx)
        if recordBreaksY(Yidx(i-1)) == 2 || Yidx(i-1) ==1
            nullclineValuesX = YnullclinesXdata(Yidx(i-1):Yidx(i)-1);
            nullclineValuesY = YnullclinesYdata(Yidx(i-1):Yidx(i)-1);

            nullclineValuesXrounded = YnullclinesXdataRounded(Yidx(i-1):Yidx(i)-1);
            nullclineValuesYrounded = YnullclinesYdataRounded(Yidx(i-1):Yidx(i)-1);
        elseif recordBreaksY(Yidx(i-1)) == 1
            nullclineValuesX = YnullclinesXdata(Yidx(i-1):Yidx(i)-1);
            nullclineValuesY = YnullclinesYdata(Yidx(i-1):Yidx(i)-1);

            nullclineValuesXrounded = YnullclinesXdataRounded(Yidx(i-1):Yidx(i)-1);
            nullclineValuesYrounded = YnullclinesYdataRounded(Yidx(i-1):Yidx(i)-1);
        end
        if length(nullclineValuesX) > shortestnullcline
            nullValues = gin(nullclineValuesX,nullclineValuesY) - nullclineValuesY;
            nr_nullValues = abs(nullValues) < varianceY;
            
            nullclineValuesX = nullclineValuesX(nr_nullValues);
            nullclineValuesY = nullclineValuesY(nr_nullValues);
            
            nullclineValuesXrounded = nullclineValuesXrounded(nr_nullValues);
            nullclineValuesYrounded = nullclineValuesYrounded(nr_nullValues);

            if length(nullclineValuesX) > shortestnullcline
                YbranchesXdata(1:length(nullclineValuesX),i-1) = nullclineValuesX;
                YbranchesYdata(1:length(nullclineValuesX),i-1) = nullclineValuesY;

                YbranchesXdataRounded(1:length(nullclineValuesX),i-1) = nullclineValuesXrounded;
                YbranchesYdataRounded(1:length(nullclineValuesX),i-1) = nullclineValuesYrounded;
            end
        end
    end

    if ~isempty(YbranchesXdata)
        idx = ~(isnan(YbranchesXdata(1,:)));
        YbranchesXdata = YbranchesXdata(:,idx);
        YbranchesYdata = YbranchesYdata(:,idx);
        
        YbranchesXdataRounded = YbranchesXdataRounded(:,idx);
        YbranchesYdataRounded = YbranchesYdataRounded(:,idx);

        nr_nullclinesY = sum(idx);
    end
   
        % combine nullclines if requested
    if options.CombineNullcline == "on"
        if options.OptimizeNullclineCombination == "on"
            if ~isempty(XbranchesXdata)
                [X_possibleCombosX, X_possibleCombosY, X_finalNrForAttempt] = AttemptNullclineCombine(XbranchesXdata, XbranchesYdata, nr_nullclinesX, max_x, min_x, max_y, min_y,1);
            end

            if ~isempty(YbranchesXdata)
                [Y_possibleCombosX, Y_possibleCombosY, Y_finalNrForAttempt] = AttemptNullclineCombine(YbranchesXdata, YbranchesYdata, nr_nullclinesY, max_x, min_x, max_y, min_y,1);
            end
        else
            if ~isempty(XbranchesXdata)
                [X_possibleCombosX, X_possibleCombosY, X_finalNrForAttempt] = AttemptNullclineCombine(XbranchesXdata, XbranchesYdata, nr_nullclinesX, max_x, min_x, max_y, min_y,0);
            end

            if ~isempty(YbranchesXdata)
                [Y_possibleCombosX, Y_possibleCombosY, Y_finalNrForAttempt] = AttemptNullclineCombine(YbranchesXdata, YbranchesYdata, nr_nullclinesY, max_x, min_x, max_y, min_y,0);
            end
        end

        if ~isempty(XbranchesXdata)
            X_best = min(X_finalNrForAttempt);
            X_bestCombo = find(X_finalNrForAttempt == X_best);
            XbranchesXdata = X_possibleCombosX{X_bestCombo(1)};
            XbranchesYdata = X_possibleCombosY{X_bestCombo(1)};
            nr_nullclinesX = X_best;
        end
        
        if ~isempty(YbranchesXdata)
            Y_best = min(Y_finalNrForAttempt);
            Y_bestCombo = find(Y_finalNrForAttempt == Y_best);
            YbranchesXdata = Y_possibleCombosX{Y_bestCombo(1)};
            YbranchesYdata = Y_possibleCombosY{Y_bestCombo(1)};
            nr_nullclinesY = Y_best;
        end
    end


    % Plot nullcline branches 
        % a) X nullclines
        % Set colours
    if nr_nullclinesX <=3
        XColours = [0.35, 0.5, 1; 0, 0, 1; 0, 0, 0.6];
    else 
        XColours = colormap(ax,winter(nr_nullclinesX+4));
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
        YColours = colormap(ax,autumn(nr_nullclinesY+4));
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
        xreal = round(xreal, 10);

        yvals = XbranchesYdata(:,i);
        yreal = yvals(~isnan(yvals));
        yreal = round(yreal,10);

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
        NIOfnctn = fimplicit(ax,NIO,[minx-PlotToleranceX maxx+PlotToleranceX miny-PlotToleranceY maxy+PlotToleranceY],'LineStyle',"none",'MeshDensity', mesh,'HandleVisibility','off');
        
        NIOxData = NIOfnctn.XData;
        NIOyData = NIOfnctn.YData;

        NIOsubValues = NIO(NIOxData, NIOyData);
        NIOtruePtsIdx = (abs(NIOsubValues) < varianceX) | isnan(NIOsubValues);
        NIOtrueDataX = NIOxData(NIOtruePtsIdx);
        NIOtrueDataY = NIOyData(NIOtruePtsIdx);

        plot(ax,NIOtrueDataX,NIOtrueDataY,'Color',colr,'LineWidth',1.5,'HandleVisibility','off');
        
        % evaluate NIO on grid
        if xfunc == 1 % function of x
            for j = i+1:nr+2:length(xaxis)
                for k = i+2:nr+2:length(yaxis)
                    symbolx = xaxis(j);
                    symboly = yaxis(k);
                    NIOvalue = NIO(symbolx,symboly);
                    if isreal(NIOvalue) && NIOvalue>0 % positive put up arrow
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','^','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    elseif isreal(NIOvalue) && NIOvalue<0 % negative put down arrow
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','v','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    else % if not defined put star
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',10,'Marker','pentagram','LineStyle','-','LineWidth',1,'HandleVisibility','off');
                    end
                end
            end
        elseif yfunc ==1 % function of y
           for j = i+1:nr+2:length(xaxis)
                for k = i+2:nr+2:length(yaxis)
                    symbolx = xaxis(j);
                    symboly = yaxis(k);
                    NIOvalue = NIO(symbolx,symboly);
                    if isreal(NIOvalue) && NIOvalue>0 % positive put right arrow
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','>','LineStyle','-','LineWidth',0.5,'HandleVisibility','off');
                    elseif isreal(NIOvalue) && NIOvalue<0 % negative put left arrow
                        plot(ax,[symbolx,symbolx],[symboly,symboly],'MarkerFaceColor',colr,'MarkerEdgeColor',colr,'MarkerSize',8,'Marker','<','LineStyle','-','LineWidth',0.5,'HandleVisibility','off');
                    else % undefined, put star
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
        xreal = round(xreal,10);

        yvals = YbranchesYdata(:,i);
        yreal = yvals(~isnan(yvals));
        yreal = round(yreal,10);

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
        NIOfnctn = fimplicit(ax,NIO,[minx-PlotToleranceX maxx+PlotToleranceX miny-PlotToleranceY maxy+PlotToleranceY],'LineStyle',"none",'MeshDensity', mesh,'HandleVisibility','off');
        
        NIOxData = NIOfnctn.XData;
        NIOyData = NIOfnctn.YData;

        NIOsubValues = NIO(NIOxData, NIOyData);
        NIOtruePtsIdx = (abs(NIOsubValues) < varianceY) | isnan(NIOsubValues);
        NIOtrueDataX = NIOxData(NIOtruePtsIdx);
        NIOtrueDataY = NIOyData(NIOtruePtsIdx);

        plot(ax,NIOtrueDataX,NIOtrueDataY,'Color',colr,'LineWidth',1.5,'HandleVisibility','off');
        
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
        fcutf = @(x,y) fin(x,y) - cutoffvalx(x,y); % subtract value from x function
        % plot when next iteration x-value is the cutoff value for x
        fimplicit(ax,fcutf,':', 'Color',[0.45 0.85 0.45], 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('f=',num2str(options.CutoffValues(1))))
        cutoffvaly = @(x,y) options.CutoffValues(2);
        gcutf = @(x,y) gin(x,y) - cutoffvaly(x,y);
        fimplicit(ax,gcutf,':', 'Color',[0 0.45 0], 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('g=',num2str(options.CutoffValues(2))))

    end

    % Lines along coordinate axes
    xline(ax,0,"LineStyle",'--',"Color",[0.5 0.5 0.5],'HandleVisibility','off')
    yline(ax,0,"LineStyle",'--',"Color",[0.5 0.5 0.5],'HandleVisibility','off')
    
    % Set axes limits
        % make a little bigger to better see nullclines along axes
    ax.XLim = [(minx-PlotToleranceX) (maxx+PlotToleranceX)];
    ax.YLim = [(miny-PlotToleranceY) (maxy+PlotToleranceY)];

    % Title, axes labels and legend
    xlabel(ax,'X');
    ylabel(ax,'Y', 'rotation', 0);
    title(ax,'Phase Portrait')

    legend(ax,"Location","northeastoutside")
    set(ax,'FontSize',18)
end

function [breaks] = splitNullclineData(XData, YData, breaks)
    diffx = diff([XData; YData],1,2); % diff between x vals in row 1 and diff between y vals in row 2
    diffx(diffx>0) = 1; % set positive differences to 1
    diffx(diffx<0) = -1; % set negative differences to -1
    
    variable = []; % if variable = 1 the nullcline is labelled as a function of x, 0 for a function of y 
    skipNext = 0; % 1 if we do not wish to examine the next difference 
    for i = 1:size(diffx,2)
        if skipNext
            skipNext = 0; % don't do anything for this difference, just switch SkipNext back to 0
        elseif isnan(diffx(1,i)) || isnan(diffx(2,i)) % NaN value
            breaks(i+1) = 3; % remove the NaN value
            breaks(i+2) = 2; % point after NaN starts a new nullcline branch
            skipNext = 1; % no need to check next difference, since it will also be NaN
            variable = [];
        elseif diffx(1,i) == 0 && diffx(2,i) == 0 %same value
            breaks(i+1) = 3; % remove value
        elseif diffx(1,i) == 0 % constant x values
            if variable == 1 % if nullcline is a function of x
                breaks(i+1) = 1; % mark the next value as the start of a new nulcline
                variable = []; % reset the variable
            else % mark the nullcline as a function of y
                variable = 0;
            end
        elseif diffx(2,i) == 0 % constant y values
            if variable == 0 % if nullcline is a function of y
                breaks(i+1) = 1; % mark the next value as the start of a new nulcline
                variable = []; % reset the variable
            else % mark the nullcline as a function of x 
                variable = 1;
            end
        elseif i>=2 && breaks(i)~=2 && breaks(i)~=1 % not the start of a new nullcline
            firstdiff_x = diffx(1,i);
            firstdiff_y = diffx(2,i);
            if breaks(i) == 3
                if breaks(i-1)==2
                    seconddiff_x = firstdiff_x;
                    seconddiff_y = firstdiff_y;
                else
                    if i > 2
                        seconddiff_x = diffx(1,i-2);
                        seconddiff_y = diffx(2,i-2);
                    else
                        seconddiff_x = firstdiff_x;
                        seconddiff_y = firstdiff_y;
                    end
                end
            else
                seconddiff_x = diffx(1,i-1);
                seconddiff_y = diffx(2,i-1);
            end
            if firstdiff_x - seconddiff_x ~= 0 && firstdiff_y - seconddiff_y ~= 0 % both x and y change direction
                breaks(i+1) = 1;
                variable = [];
            elseif firstdiff_x - seconddiff_x ~= 0 % x values change direction
                if variable == 1
                    breaks(i+1) = 1;
                    variable = [];
                else
                    variable = 0;
                end
            elseif firstdiff_y - seconddiff_y ~= 0 % y values change direction
                if variable == 0
                    breaks(i+1) = 1;
                    variable = [];
                else
                    variable = 1;
                end
            end
        end
    end
end

function [combinationOptionsX, combinationOptionsY, nrNullclinesOptions] = AttemptNullclineCombine(branchDataX, branchDataY, nrNullclines, XMax, XMin, YMax, YMin, optimize)
        nullclineFirstLast = zeros(nrNullclines*2,2);

        % stores a vector form of the direction at each end point
        vectorFormsofEnds = zeros(nrNullclines*2,2);
        
        % first values
        firstX = branchDataX(1,:); 
        firstX = firstX.';
        nullclineFirstLast(1:nrNullclines,1) = firstX;
        firstY = branchDataY(1,:);
        firstY = firstY.';
        nullclineFirstLast(1:nrNullclines,2) = firstY;

        for i = 1:nrNullclines
            % end values
            branchX = branchDataX(1:end,i); % extract X values for nullcline branch
            branchX = branchX(~isnan(branchX)); % remove nan values
            branchXend = branchX(end); % take second last X value which is a number
            nullclineFirstLast(nrNullclines+i,1) = branchXend; % save to array
            branchY = branchDataY(1:end,i);
            branchY = branchY(~isnan(branchY));
            branchYend = branchY(end);
            nullclineFirstLast(nrNullclines+i,2) = branchYend;

            % vector forms of ends of nullclines
            buffer = floor(length(branchX)/10); 
                    % how far into the branch to consider when taking the direction vector at the end
                    % with error taking adjacent points may not be accurate
                    % so we average over a larger section

            branchX_secondLast = branchX(end-buffer);
            branchY_secondLast = branchY(end-buffer);

            % take difference between end and the value slightly away from
            % the end
            diffLast2_X = branchXend - branchX_secondLast;
            diffLast2_Y = branchYend - branchY_secondLast;
            vectorFormsofEnds(i+nrNullclines,1) = diffLast2_X;
            vectorFormsofEnds(i+nrNullclines,2) = diffLast2_Y;

            % vector forms of starts of nullclines
            branchX_second = branchX(1+buffer);
            branchY_second = branchY(1+buffer);

            diffFirst2_X = branchX_second - nullclineFirstLast(i,1);
            diffFirst2_Y = branchY_second - nullclineFirstLast(i,2);

            vectorFormsofEnds(i,1) = diffFirst2_X;
            vectorFormsofEnds(i,2) = diffFirst2_Y;
        end

        % vector magnitudes
        squared = vectorFormsofEnds.^2;
        summed  = squared(:,1) + squared(:,2);
        vectorMagnitudes = sqrt(summed);

        % find distances between points 
        distances = pdist(nullclineFirstLast);
        distances = squareform(distances); % puts distances in square matrix

        % min distance to say two endpoints are close enough to try
        % combining -- can be changed to make the method more or less
        % specific
        closenessValue = 10;
        distmin = sqrt(((XMax-XMin)/(closenessValue))^2 + ((YMax-YMin)/(closenessValue))^2);

        % set distances we don't want to evaluate to greater than the min
        % distance 
            % we don't need to compare the distance of a point to itself
            % if we test the distance from 1 to 2, we dont need the
            % distance from 2 to 1
        for i = 1:2*nrNullclines
            for j = 1:2*nrNullclines
                if i >= j || j == i + nrNullclines
                    distances(i,j) = distmin + 1;
                end
            end
        end

        %find distances smaller than min distance 
        [iIdx, jIdx] = find(distances <= distmin);
        comparisons = [iIdx jIdx];

        if optimize == 1
            if size(comparisons,1) <= 10
                nrs = 1:size(comparisons,1);
                combinationOrders = perms(nrs);
            else
                combinationOrders = zeros(size(comparisons,1),size(comparisons,1));
                nrs = 1:size(comparisons,1);
                for s = 1:size(comparisons,1)
                    combinationOrders(s,:) = circshift(nrs,s);
                end
            end
        else
            combinationOrders = 1:size(comparisons,1);
        end

        combinationAttempts = size(combinationOrders,1);
        
        nrNullclinesOptions = NaN(combinationAttempts,1);
        combinationOptionsX = cell(combinationAttempts,1);
        combinationOptionsY = cell(combinationAttempts,1);

        for g = 1:combinationAttempts
            comparisons2 = comparisons(combinationOrders(g,:),:);

            branchDataXduplicate = branchDataX;
            branchDataYduplicate = branchDataY;

            nullclinesToDelete = [];

            while ~isempty(comparisons2)
                
                % we start with a point and find all the points near it
                pair1 = comparisons2(1,1);
                [AllNear, ~] = find(comparisons2 == pair1);
                Allcomparisons = comparisons2(AllNear,:);
                AllNear = [Allcomparisons(:,1); Allcomparisons(:,2)];
                AllNear = unique(AllNear);
                
                % all the pairings in this local area 
                pairings = nchoosek(AllNear, 2);
    
                % initialize for storing angles between pairs of end points 
                anglesNear = zeros(size(pairings,1),1);
    
                for i = 1:size(pairings,1)
                    val1 = pairings(i,1); % first endpoint
                    val2 = pairings(i,2); % second endpoint
                    vector1 = vectorFormsofEnds(val1,:);
                    vector2 = vectorFormsofEnds(val2,:);
                    
                    % take dot product of vector forms of the ends divided by
                    % the magnitudes of the vectors, then take the inverse cos
                    % to get the angle between the ends
                    % vectors are all pointing away from the end, so, if the
                    % direction of the ends is aligned it will have a value
                    % close to pi
                    vectorDot = dot(vector1,vector2);
                    angle = vectorDot/(vectorMagnitudes(val1)*vectorMagnitudes(val2));
                    angle = acos(angle);
                    anglesNear(i) = angle;
                end
    
                totalPairings = floor(length(AllNear)/2); % ex with 3 points we can make 1 pairing
                while totalPairings > 0 && ~isempty(pairings)
                    maxAngle = max(anglesNear);
                    bestPairIdx = find(anglesNear == maxAngle); 
                    bestPair = pairings(bestPairIdx(1),:); % find the pair of endpoints where the angle between the ends is the closest to pi
                    point1 = bestPair(1);
                    point2 = bestPair(2);
                    
                    % extract data 
                    if point1 > nrNullclines % indicates point is at the end of the nullcline in column point1 - nr_nullclinesX
                        branch1X = branchDataXduplicate(:,point1 - nrNullclines);
                        branch1X = branch1X(~isnan(branch1X));
                        branch1Y = branchDataYduplicate(:,point1 - nrNullclines);
                        branch1Y = branch1Y(~isnan(branch1Y));
                    else 
                        branch1X = branchDataXduplicate(:,point1);
                        branch1X = branch1X(~isnan(branch1X));
                        branch1Y = branchDataYduplicate(:,point1);
                        branch1Y = branch1Y(~isnan(branch1Y));
                    end
    
                    if point2 > nrNullclines
                        branch2X = branchDataXduplicate(:,point2 - nrNullclines);
                        branch2X = branch2X(~isnan(branch2X));
                        branch2Y = branchDataYduplicate(:,point2 - nrNullclines);
                        branch2Y = branch2Y(~isnan(branch2Y));
                    else
                        branch2X = branchDataXduplicate(:,point2);
                        branch2X = branch2X(~isnan(branch2X));
                        branch2Y = branchDataYduplicate(:,point2);
                        branch2Y = branch2Y(~isnan(branch2Y));
                    end
    
                    % try to combine
                    flip1 = 0; % 1 if we flip the direction of the first nullcline
                    flip2 = 0; % 1 if we flip the direction of the second nullcline
                    
                    if point1 <= nrNullclines %if point 1 is a beginning point
                        % flip so that the beginning point is at the end, to be
                        % combined with the second nullcline
                        branch1X = flip(branch1X);
                        branch1Y = flip(branch1Y);
                        flip1 = 1;
                    end
    
                    if point2 > nrNullclines % if point 2 is the end point of its nullcline
                        % flip numbers so that the end point is at the
                        % beginning next to the values for the first nullcline
                        branch2X = flip(branch2X);
                        branch2Y = flip(branch2Y);
                        flip2 = 1;                   
                    end
    
                    % concatenate values
                    comboAttemptX = [branch1X; branch2X];
                    comboAttemptY = [branch1Y; branch2Y];
                    
                    % take difference between points in the new combination
                    xDifference = diff(comboAttemptX);
                    yDifference = diff(comboAttemptY);
    
                    testX = xDifference == 0;
                    testY = yDifference == 0;
                    duplicate = testX & testY;
                    xDifference(duplicate) = [];
                    yDifference(duplicate) = [];
    
                    comboAttemptX(duplicate) = [];
                    comboAttemptY(duplicate) = [];
    
    
                    % if it is consistently increasing or decreasing in either
                    % one of the variables, then it is a valid combination
                    if all(xDifference>0) || all(xDifference<0) || all(yDifference>0) || all(yDifference<0)
                        % add combined nullcline to the column associated with
                        % the nullcline for point 1
                        if point1 <= nrNullclines % if point1 a start value
                            branchDataXduplicate(1:length(comboAttemptX),point1) = comboAttemptX;
                            branchDataYduplicate(1:length(comboAttemptY),point1) = comboAttemptY; %reflect changes in rounded version?
                        else
                            branchDataXduplicate(1:length(comboAttemptX),point1-nrNullclines) = comboAttemptX.';
                            branchDataYduplicate(1:length(comboAttemptY),point1-nrNullclines) = comboAttemptY.';
                        end
                        
                        % add second nullcline to delete list
                        if point2 <= nrNullclines
                             nullclinesToDelete = [nullclinesToDelete; point2];
                        else
                             nullclinesToDelete = [nullclinesToDelete; point2-nrNullclines];
                        end
    
    
                        % delete instances of either point in lists of local
                        % pairings and total list of comparisons
                        totalPairings = totalPairings - 1;
                        todeletePairs = any((pairings == point1)|(pairings == point2),2);
                        pairings(todeletePairs,:) = [];
    
                        anglesNear(todeletePairs)=[];
    
                        todeleteComparisons = any((comparisons2 == point1)|(comparisons2 == point2),2);
                        comparisons2(todeleteComparisons,:) = [];
    
                        % if there are still comparisons left, we swap the
                        % labels on the endpoints to reflect the changes due to
                        % the nullcline combination
                        if ~isempty(comparisons2) 
                            if flip1==1 
                                vectorFormsofEnds(nrNullclines+point1,:) = vectorFormsofEnds(point1,:);
    
                                % when in pairings it says end of 1 put start of 1
                                switchlabels = (pairings == nrNullclines + point1);
                                pairings(switchlabels) = point1;
    
                                switchlabels2 = (comparisons2 == nrNullclines + point1);
                                comparisons2(switchlabels2) = point1;
                            end
                            if flip2==1
                                % when it says start of 2 replace with end of 1
                                switchlabels = (pairings == point2 - nrNullclines);
                                switchlabels2 = (comparisons2 == point2 - nrNullclines);
                                if flip1==1
                                    vectorFormsofEnds(point2 - nrNullclines,:) = vectorFormsofEnds(point1 + nrNullclines,:);
    
                                    pairings(switchlabels) = point1  + nrNullclines;
                                    comparisons2(switchlabels2) = point1  + nrNullclines;
                                else
                                    vectorFormsofEnds(point2 - nrNullclines,:) = vectorFormsofEnds(point1,:);
    
                                    pairings(switchlabels) = point1;
                                    comparisons2(switchlabels2) = point1;
                                end
                            else 
                                % when it says end of 2 replace with end of 1
                                switchlabels = (pairings == point2 + nrNullclines);
                                switchlabels2 = (comparisons2 == point2 + nrNullclines);
                                if flip1==1
                                    vectorFormsofEnds(point2 + nrNullclines,:) = vectorFormsofEnds(point1 + nrNullclines,:);
    
                                    pairings(switchlabels) = point1 + nrNullclines;
                                    comparisons2(switchlabels2) = point1 + nrNullclines;
                                else
                                    vectorFormsofEnds(point2 + nrNullclines,:) = vectorFormsofEnds(point1,:);
    
                                    pairings(switchlabels) = point1;
                                    comparisons2(switchlabels2) = point1;
                                end
    
                            end
                        end                                                
                    else
                        % remove this combo from pairs and comparisons
                        pairings(bestPairIdx(1),:) = [];
                        
                        row1 = [point1 point2];
                        row2 = [point2 point1];
    
                        idx1 = ismember(comparisons2,row1,"row");
                        idx2 = ismember(comparisons2,row2,"row");
    
                        % if both points appear in the combination in either
                        % order, we remove this comparison
                        comparisons2(idx1 | idx2,:) = [];
    
                        % remove angle between this comparison
                        anglesNear(bestPairIdx(1)) = [];
    
                    end
                end
    
            end
            % delete extra nullclines
            branchDataXduplicate(:,nullclinesToDelete) = [];
            branchDataYduplicate(:,nullclinesToDelete) = [];
    
            % update number of nullclines
            nrNullclinesOptions(g) = size(branchDataXduplicate,2);

            combinationOptionsX{g} = branchDataXduplicate;
            combinationOptionsY{g} = branchDataYduplicate;

        end
end

    