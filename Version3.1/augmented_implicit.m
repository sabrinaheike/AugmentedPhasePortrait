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
    PlotToleranceX = (max_x-min_x)/1000; % fimplicit looks for nullclines in a slightly larger range than the one inputted so that nullclines on the edges of the range are not missed. 
    PlotToleranceY = (max_y-min_y)/1000;
    XnullclinesPlot = fimplicit(ax,Xnullclines,[min_x-PlotToleranceX max_x+PlotToleranceX min_y-PlotToleranceY max_y+PlotToleranceY], 'LineStyle', 'none', 'HandleVisibility','off','MeshDensity',300);
    hold(ax,"on")
    YnullclinesPlot = fimplicit(ax,Ynullclines,[min_x-PlotToleranceX max_x+PlotToleranceX min_y-PlotToleranceY max_y+PlotToleranceY], 'LineStyle', 'none', 'HandleVisibility','off','MeshDensity',300);

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
    skipNext = 0; % 1 if we do not wish to examine the next difference 
    for i = 1:size(diffx,2)
        if skipNext
            skipNext = 0; % don't do anything for this difference, just switch SkipNext back to 0
        elseif isnan(diffx(1,i)) || isnan(diffx(2,i)) % NaN value
            recordBreaksX(i+1) = 3; % remove the NaN value
            recordBreaksX(i+2) = 2; % point after NaN starts a new nullcline branch
            skipNext = 1; % no need to check next difference, since it will also be NaN
            variable = [];
        elseif diffx(1,i) == 0 && diffx(2,i) == 0 %same value
            recordBreaksX(i+1) = 3; % remove value
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
        elseif i>=2 && recordBreaksX(i)~=2 % not the start of a new nullcline
            firstdiff_x = diffx(1,i);
            firstdiff_y = diffx(2,i);
            if recordBreaksX(i) == 3
                if recordBreaksX(i-1)==2
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
    XnullclinesXdata = XnullclinesXdata(recordBreaksX~=3);
    XnullclinesYdata = XnullclinesYdata(recordBreaksX~=3);

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

        % sort the data into the array 
    shortestnullcline = 20; % don't want to do calculations for nullcline segments which are too short to properly see
    variance = (maxx-minx)/100; % how far we allow the x-values to change in the next iteration from a nullcline point
    for i = 2:length(Xidx)
        if recordBreaksX(Xidx(i-1)) == 2 || Xidx(i-1) == 1 % don't share a common point between nullcline branches
            nullclineValuesX = XnullclinesXdata(Xidx(i-1):Xidx(i)-1);
            nullclineValuesY = XnullclinesYdata(Xidx(i-1):Xidx(i)-1);
        elseif recordBreaksX(Xidx(i-1)) == 1 % share a common point between nullcline branches, new branch starts due to change in direction not becuase of a NaN
            nullclineValuesX = XnullclinesXdata(Xidx(i-1)-1:Xidx(i)-1);
            nullclineValuesY = XnullclinesYdata(Xidx(i-1)-1:Xidx(i)-1);
        end
        if length(nullclineValuesX) > shortestnullcline % if the nullcline is long enough
            
            % for each point on the nullcline, take the absolute difference
            % between the x value of the point, and the x value of the
            % point that it is mapped to. If it is a true nullcline this
            % should be very close to 0 (within our threshold of variance)
            nullValues = fin(nullclineValuesX,nullclineValuesY) - nullclineValuesX; 
            nr_nullValues = abs(nullValues) < variance; % indices of true nullcline values
            nullclineValuesX = nullclineValuesX(nr_nullValues);
            nullclineValuesY = nullclineValuesY(nr_nullValues);
            if length(nullclineValuesX) > shortestnullcline %if the nullcline with erroneous points removed is still long enough
                % add values to array 
                XbranchesXdata(1:length(nullclineValuesX),i-1) = nullclineValuesX;
                XbranchesYdata(1:length(nullclineValuesX),i-1) = nullclineValuesY;
            end
        end
    end

    % update the number of X-nullclines since it might have changed when we
    % removed nullclines that were too short or were not true nullclines
    idx = ~(isnan(XbranchesXdata(1,:)));
    XbranchesXdata = XbranchesXdata(:,idx);
    XbranchesYdata = XbranchesYdata(:,idx);
    nr_nullclinesX = sum(idx);

        % b) Y nullclines

    diffy = diff([YnullclinesXdata; YnullclinesYdata],1,2);
    diffy(diffy>0) = 1;
    diffy(diffy<0) = -1;
    recordBreaksY = zeros(length(YnullclinesXdata),1);
    variable = [];
    skipNext = 0;
    for i = 1:size(diffy,2)
        if skipNext
            skipNext = 0;
        elseif isnan(diffy(1,i)) || isnan(diffy(2,i))
            recordBreaksY(i+1) = 3;
            recordBreaksY(i+2) = 2;
            skipNext = 1;
            variable = [];
        elseif diffy(1,i) == 0 && diffy(2,i) == 0 %same value
            recordBreaksY(i+1) = 3; % remove value
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
        elseif i>=2 && recordBreaksY(i)~=2
            firstdiff_x = diffy(1,i);
            firstdiff_y = diffy(2,i);
            if recordBreaksY(i) == 3
                if recordBreaksY(i-1)==2
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
    YnullclinesXdata = YnullclinesXdata(recordBreaksY~=3);
    YnullclinesYdata = YnullclinesYdata(recordBreaksY~=3);

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

    YbranchesXdata = nan(length(YnullclinesXdata),nr_nullclinesY);
    YbranchesYdata = nan(length(YnullclinesXdata),nr_nullclinesY);

    % sort the data into the array
    variance = (maxy-miny)/100;
    for i = 2:length(Yidx)
        if recordBreaksY(Yidx(i-1)) == 2 || Yidx(i-1) ==1
            nullclineValuesX = YnullclinesXdata(Yidx(i-1):Yidx(i)-1);
            nullclineValuesY = YnullclinesYdata(Yidx(i-1):Yidx(i)-1);
        elseif recordBreaksY(Yidx(i-1)) == 1
            nullclineValuesX = YnullclinesXdata(Yidx(i-1)-1:Yidx(i)-1);
            nullclineValuesY = YnullclinesYdata(Yidx(i-1)-1:Yidx(i)-1);
        end
        if length(nullclineValuesX) > shortestnullcline
            nullValues = gin(nullclineValuesX,nullclineValuesY) - nullclineValuesY;
            nr_nullValues = abs(nullValues) < variance;
            nullclineValuesX = nullclineValuesX(nr_nullValues);
            nullclineValuesY = nullclineValuesY(nr_nullValues);
            if length(nullclineValuesX) > shortestnullcline
                YbranchesXdata(1:length(nullclineValuesX),i-1) = nullclineValuesX;
                YbranchesYdata(1:length(nullclineValuesX),i-1) = nullclineValuesY;
            end
        end
    end

    idx = ~(isnan(YbranchesXdata(1,:)));
    YbranchesXdata = YbranchesXdata(:,idx);
    YbranchesYdata = YbranchesYdata(:,idx);
    nr_nullclinesY = sum(idx);
   
        % combine nullclines if requested
    if options.CombineNullcline == "on"
        % stores both end points of each nullcline
        nullclineFirstLastX = zeros(nr_nullclinesX*2,2);

        % stores a vector form of the direction at each end point
        vectorFormsofEndsX = zeros(nr_nullclinesX*2,2);
        
        % first values
        XfirstX = XbranchesXdata(2,:); % take second point because sometimes there is a numerical error right near the end
        XfirstX = XfirstX.';
        nullclineFirstLastX(1:nr_nullclinesX,1) = XfirstX;
        XfirstY = XbranchesYdata(2,:);
        XfirstY = XfirstY.';
        nullclineFirstLastX(1:nr_nullclinesX,2) = XfirstY;

        for i = 1:nr_nullclinesX
            % end values
            XbranchX = XbranchesXdata(1:end,i); % extract X values for nullcline branch
            XbranchX = XbranchX(~isnan(XbranchX)); % remove nan values
            XbranchXend = XbranchX(end-1); % take second last X value which is a number
            nullclineFirstLastX(nr_nullclinesX+i,1) = XbranchXend; % save to array
            XbranchY = XbranchesYdata(1:end,i);
            XbranchY = XbranchY(~isnan(XbranchY));
            XbranchYend = XbranchY(end-1);
            nullclineFirstLastX(nr_nullclinesX+i,2) = XbranchYend;

            % vector forms of ends of nullclines
            bufferX = floor(length(XbranchX)/10); 
                    % how far into the branch to consider when taking the direction vector at the end
                    % with error taking adjacent points may not be accurate
                    % so we average over a larger section

            XbranchX_secondLast = XbranchX(end-bufferX);
            XbranchY_secondLast = XbranchY(end-bufferX);

            % take difference between end and the value slightly away from
            % the end
            XdiffLast2_X = XbranchXend - XbranchX_secondLast;
            XdiffLast2_Y = XbranchYend - XbranchY_secondLast;
            vectorFormsofEndsX(i+nr_nullclinesX,1) = XdiffLast2_X;
            vectorFormsofEndsX(i+nr_nullclinesX,2) = XdiffLast2_Y;

            % vector forms of starts of nullclines
            XbranchX_second = XbranchX(1+bufferX);
            XbranchY_second = XbranchY(1+bufferX);

            XdiffFirst2_X = XbranchX_second - nullclineFirstLastX(i,1);
            XdiffFirst2_Y = XbranchY_second - nullclineFirstLastX(i,2);

            vectorFormsofEndsX(i,1) = XdiffFirst2_X;
            vectorFormsofEndsX(i,2) = XdiffFirst2_Y;
        end

        % vector magnitudes
        squared = vectorFormsofEndsX.^2;
        summed  = squared(:,1) + squared(:,2);
        vectorMagnitudesX = sqrt(summed);

        % find distances between points 
        distancesX = pdist(nullclineFirstLastX);
        distancesX = squareform(distancesX); % puts distances in square matrix

        % min distance to say two endpoints are close enough to try
        % combining -- can be changed to make the method more or less
        % specific
        closenessValue = 10;
        distmin = sqrt(((max_x-min_x)/(closenessValue))^2 + ((max_y-min_y)/(closenessValue))^2);

        % set distances we don't want to evaluate to greater than the min
        % distance 
            % we don't need to compare the distance of a point to itself
            % if we test the distance from 1 to 2, we dont need the
            % distance from 2 to 1
        for i = 1:2*nr_nullclinesX
            for j = 1:2*nr_nullclinesX
                if i >= j || j == i + nr_nullclinesX
                    distancesX(i,j) = distmin + 1;
                end
            end
        end

        %find distances smaller than min distance 
        [iIdx, jIdx] = find(distancesX <= distmin);
        comparisonsX = [iIdx jIdx];
        nullclinesToDeleteX = [];

        while ~isempty(comparisonsX)
            
            % we start with a point and find all the points near it
            pair1 = comparisonsX(1,1);
            [AllNearX, ~] = find(comparisonsX == pair1);
            Allcomparisons = comparisonsX(AllNearX,:);
            AllNear = [Allcomparisons(:,1); Allcomparisons(:,2)];
            AllNear = unique(AllNear);
            
            % all the pairings in this local area 
            pairings = nchoosek(AllNear, 2);

            % initialize for storing angles between pairs of end points 
            anglesNear = zeros(size(pairings,1),1);

            for i = 1:size(pairings,1)
                val1 = pairings(i,1); % first endpoint
                val2 = pairings(i,2); % second endpoint
                vector1 = vectorFormsofEndsX(val1,:);
                vector2 = vectorFormsofEndsX(val2,:);
                
                % take dot product of vector forms of the ends divided by
                % the magnitudes of the vectors, then take the inverse cos
                % to get the angle between the ends
                % vectors are all pointing away from the end, so, if the
                % direction of the ends is aligned it will have a value
                % close to pi
                vectorDot = dot(vector1,vector2);
                angle = vectorDot/(vectorMagnitudesX(val1)*vectorMagnitudesX(val2));
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
                if point1 > nr_nullclinesX % indicates point is at the end of the nullcline in column point1 - nr_nullclinesX
                    branch1X = XbranchesXdata(:,point1 - nr_nullclinesX);
                    branch1X = branch1X(~isnan(branch1X));
                    branch1Y = XbranchesYdata(:,point1 - nr_nullclinesX);
                    branch1Y = branch1Y(~isnan(branch1Y));
                else 
                    branch1X = XbranchesXdata(:,point1);
                    branch1X = branch1X(~isnan(branch1X));
                    branch1Y = XbranchesYdata(:,point1);
                    branch1Y = branch1Y(~isnan(branch1Y));
                end

                if point2 > nr_nullclinesX
                    branch2X = XbranchesXdata(:,point2 - nr_nullclinesX);
                    branch2X = branch2X(~isnan(branch2X));
                    branch2Y = XbranchesYdata(:,point2 - nr_nullclinesX);
                    branch2Y = branch2Y(~isnan(branch2Y));
                else
                    branch2X = XbranchesXdata(:,point2);
                    branch2X = branch2X(~isnan(branch2X));
                    branch2Y = XbranchesYdata(:,point2);
                    branch2Y = branch2Y(~isnan(branch2Y));
                end

                % try to combine
                flip1 = 0; % 1 if we flip the direction of the first nullcline
                flip2 = 0; % 1 if we flip the direction of the second nullcline
                
                if point1 <= nr_nullclinesX %if point 1 is a beginning point
                    % flip so that the beginning point is at the end, to be
                    % combined with the second nullcline
                    branch1X = flip(branch1X);
                    branch1Y = flip(branch1Y);
                    flip1 = 1;
                end

                if point2 > nr_nullclinesX % if point 2 is the end point of its nullcline
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

                % if it is consistently increasing or decreasing in either
                % one of the variables, then it is a valid combination
                if all(xDifference>0) || all(xDifference<0) || all(yDifference>0) || all(yDifference<0)
                    % add combined nullcline to the column associated with
                    % the nullcline for point 1
                    if point1 <= nr_nullclinesX % if point1 a start value
                        XbranchesXdata(1:length(comboAttemptX),point1) = comboAttemptX;
                        XbranchesYdata(1:length(comboAttemptY),point1) = comboAttemptY;
                    else
                        XbranchesXdata(1:length(comboAttemptX),point1-nr_nullclinesX) = comboAttemptX.';
                        XbranchesYdata(1:length(comboAttemptY),point1-nr_nullclinesX) = comboAttemptY.';
                    end
                    
                    % add second nullcline to delete list
                    if point2 <= nr_nullclinesX
                         nullclinesToDeleteX = [nullclinesToDeleteX; point2];
                    else
                         nullclinesToDeleteX = [nullclinesToDeleteX; point2-nr_nullclinesX];
                    end


                    % delete instances of either point in lists of local
                    % pairings and total list of comparisons
                    totalPairings = totalPairings - 1;
                    todeletePairs = any((pairings == point1)|(pairings == point2),2);
                    pairings(todeletePairs,:) = [];

                    anglesNear(todeletePairs)=[];

                    todeleteComparisons = any((comparisonsX == point1)|(comparisonsX == point2),2);
                    comparisonsX(todeleteComparisons,:) = [];

                    % if there are still comparisons left, we swap the
                    % labels on the endpoints to reflect the changes due to
                    % the nullcline combination
                    if ~isempty(comparisonsX) 
                        if flip1==1 
                            % when in pairings it says end of 1 put start of 1
                            switchlabels = (pairings == nr_nullclinesX + point1);
                            pairings(switchlabels) = point1;

                            switchlabels2 = (comparisonsX == nr_nullclinesX + point1);
                            comparisonsX(switchlabels2) = point1;
                        end
                        if flip2==1
                            % when it says start of 2 replace with end of 1
                            switchlabels = (pairings == nr_nullclinesX - point2);
                            switchlabels2 = (comparisonsX == nr_nullclinesX - point2);
                            if flip1==1
                                pairings(switchlabels) = point1  + nr_nullclinesX;
                                comparisonsX(switchlabels2) = point1  + nr_nullclinesX;
                            else
                                pairings(switchlabels) = point1;
                                comparisonsX(switchlabels2) = point1;
                            end
                        else 
                            % when it says end of 2 replace with end of 1
                            switchlabels = (pairings == point2 + nr_nullclinesX);
                            switchlabels2 = (comparisonsX == point2 + nr_nullclinesX);
                            if flip1==1
                                pairings(switchlabels) = point1 + nr_nullclinesX;
                                comparisonsX(switchlabels2) = point1 + nr_nullclinesX;
                            else
                                pairings(switchlabels) = point1;
                                comparisonsX(switchlabels2) = point1;
                            end

                        end
                    end                                                
                else
                    % remove this combo from pairs and comparisons
                    pairings(bestPairIdx(1),:) = [];
                    
                    row1 = [point1 point2];
                    row2 = [point2 point1];

                    idx1 = ismember(comparisonsX,row1,"row");
                    idx2 = ismember(comparisonsX,row2,"row");

                    % if both points appear in the combination in either
                    % order, we remove this comparison
                    comparisonsX(idx1 | idx2,:) = [];

                    % remove angle between this comparison
                    anglesNear(bestPairIdx(1)) = [];

                end
            end

        end
        % delete extra nullclines
        XbranchesXdata(:,nullclinesToDeleteX) = [];
        XbranchesYdata(:,nullclinesToDeleteX) = [];

        % update number of nullclines
        nr_nullclinesX = size(XbranchesXdata,2);
    
        %%% Y nullclines (same process for Y as for X nullclines)
        nullclineFirstLastY = zeros(nr_nullclinesY*2,2);
        vectorFormsofEndsY = zeros(nr_nullclinesY*2,2);

        % first values 
        YfirstX = YbranchesXdata(2,:);
        YfirstX = YfirstX.';
        nullclineFirstLastY(1:nr_nullclinesY,1) = YfirstX;
        YfirstY = YbranchesYdata(2,:);
        YfirstY = YfirstY.';
        nullclineFirstLastY(1:nr_nullclinesY,2) = YfirstY;

        for j = 1:nr_nullclinesY
            % end values 
            YbranchX = YbranchesXdata(1:end,j);
            YbranchX = YbranchX(~isnan(YbranchX));
            YbranchXend = YbranchX(end-1);
            nullclineFirstLastY(nr_nullclinesY+j,1) = YbranchXend;

            YbranchY = YbranchesYdata(1:end,j);
            YbranchY = YbranchY(~isnan(YbranchY));
            YbranchYend = YbranchY(end-1);
            nullclineFirstLastY(nr_nullclinesY+j,2) = YbranchYend;

            % vector forms of ends of nullclines
            bufferY = floor(length(YbranchX)/10);
            
            YbranchX_secondLast = YbranchX(end-bufferY);
            YbranchY_secondLast = YbranchY(end-bufferY);

            YdiffLast2_X = YbranchXend - YbranchX_secondLast;
            YdiffLast2_Y = YbranchYend - YbranchY_secondLast;
            vectorFormsofEndsY(j+nr_nullclinesY,1) = YdiffLast2_X;
            vectorFormsofEndsY(j+nr_nullclinesY,2) = YdiffLast2_Y;

            % vector forms of starts of nullclines
            YbranchX_second = YbranchX(1+bufferY);
            YbranchY_second = YbranchY(1+bufferY);

            YdiffFirst2_X = YbranchX_second - nullclineFirstLastY(j,1);
            YdiffFirst2_Y = YbranchY_second - nullclineFirstLastY(j,2);

            vectorFormsofEndsY(j,1) = YdiffFirst2_X;
            vectorFormsofEndsY(j,2) = YdiffFirst2_Y;
        end

        % vector magnitudes
        squared = vectorFormsofEndsY.^2;
        summed = squared(:,1) + squared(:,2);
        vectorMagnitudesY = sqrt(summed);

        % find distances between points 
        distancesY = pdist(nullclineFirstLastY);
        distancesY = squareform(distancesY);

        % set distances we don't want to consider to be greater than the
        % min distance
        for i = 1:2*nr_nullclinesY
            for j = 1:2*nr_nullclinesY
                if i >= j || j == i + nr_nullclinesY
                    distancesY(i,j) = distmin + 1;
                end
            end
        end

        % find distances smaller than min distance
        [iIdx, jIdx] = find(distancesY <= distmin);
        comparisonsY = [iIdx jIdx];
        nullclinesToDeleteY = [];
        while ~isempty(comparisonsY)
            pair1 = comparisonsY(1,1);
            [AllNearY, ~] = find(comparisonsY == pair1);
            Allcomparisons = comparisonsY(AllNearY,:);
            AllNear = [Allcomparisons(:,1); Allcomparisons(:,2)];
            AllNear = unique(AllNear);
            
            pairingsY = nchoosek(AllNear, 2);

            anglesNear = zeros(size(pairingsY,1),1);

            for i = 1:size(pairingsY,1)
                val1 = pairingsY(i,1);
                val2 = pairingsY(i,2);
                vector1 = vectorFormsofEndsY(val1,:);
                vector2 = vectorFormsofEndsY(val2,:);
                vectorDot = dot(vector1,vector2);
                angle = vectorDot/(vectorMagnitudesY(val1)*vectorMagnitudesY(val2));
                angle = acos(angle);
                anglesNear(i) = angle;
            end

            totalPairings = floor(length(AllNear)/2);
            while totalPairings > 0 && ~isempty(pairingsY)
                maxAngle = max(anglesNear);
                bestPairIdx = find(anglesNear == maxAngle);
                bestPair = pairingsY(bestPairIdx(1),:);
                point1 = bestPair(1);
                point2 = bestPair(2);
                
                % extract data 
                if point1 > nr_nullclinesY
                    branch1X = YbranchesXdata(:,point1 - nr_nullclinesY);
                    branch1X = branch1X(~isnan(branch1X));
                    branch1Y = YbranchesYdata(:,point1 - nr_nullclinesY);
                    branch1Y = branch1Y(~isnan(branch1Y));
                else
                    branch1X = YbranchesXdata(:,point1);
                    branch1X = branch1X(~isnan(branch1X));
                    branch1Y = YbranchesYdata(:,point1);
                    branch1Y = branch1Y(~isnan(branch1Y));
                end

                if point2 > nr_nullclinesY
                    branch2X = YbranchesXdata(:,point2 - nr_nullclinesY);
                    branch2X = branch2X(~isnan(branch2X));
                    branch2Y = YbranchesYdata(:,point2 - nr_nullclinesY);
                    branch2Y = branch2Y(~isnan(branch2Y));
                else
                    branch2X = YbranchesXdata(:,point2);
                    branch2X = branch2X(~isnan(branch2X));
                    branch2Y = YbranchesYdata(:,point2);
                    branch2Y = branch2Y(~isnan(branch2Y));
                end

                % try to combine
                flip1 = 0;
                flip2 = 0;
                
                if point1 <= nr_nullclinesY
                    branch1X = flip(branch1X);
                    branch1Y = flip(branch1Y);
                    flip1 = 1;
                end

                if point2 > nr_nullclinesY
                    branch2X = flip(branch2X);
                    branch2Y = flip(branch2Y);
                    flip2 = 1;                   
                end

                comboAttemptX = [branch1X; branch2X];
                comboAttemptY = [branch1Y; branch2Y];
                
                xDifference = diff(comboAttemptX);
                yDifference = diff(comboAttemptY);

                if all(xDifference>0) || all(xDifference<0) || all(yDifference>0) || all(yDifference<0)
                    if point1 <= nr_nullclinesX
                        YbranchesXdata(1:length(comboAttemptX),point1) = comboAttemptX;
                        YbranchesYdata(1:length(comboAttemptY),point1) = comboAttemptY;
                    else
                        YbranchesXdata(1:length(comboAttemptX),point1-nr_nullclinesY) = comboAttemptX.';
                        YbranchesYdata(1:length(comboAttemptY),point1-nr_nullclinesY) = comboAttemptY.';
                    end
                    
                    % add second nullcline to delete list
                    if point2 <= nr_nullclinesY
                         nullclinesToDeleteY = [nullclinesToDeleteY; point2];
                    else
                         nullclinesToDeleteY = [nullclinesToDeleteY; point2-nr_nullclinesY];
                    end


                    % delete instances of either point in lists of pairings
                    totalPairings = totalPairings - 1;
                    todeletePairs = any((pairingsY == point1)|(pairingsY == point2),2);
                    pairingsY(todeletePairs,:) = [];

                    anglesNear(todeletePairs)=[];

                    todeleteComparisons = any((comparisonsY == point1)|(comparisonsY == point2),2);
                    comparisonsY(todeleteComparisons,:) = [];

                    if ~isempty(comparisonsY) % also swap comparisons? 
                        if flip1==1 
                            % when in pairings it says end of 1 put start of 1
                            switchlabels = (pairingsY == nr_nullclinesY + point1);
                            pairingsY(switchlabels) = point1;

                            switchlabels2 = (comparisonsY == nr_nullclinesY + point1);
                            comparisonsY(switchlabels2) = point1;
                        end
                        if flip2==1
                            % when it says start of 2 replace with end of 1
                            switchlabels = (pairingsY == nr_nullclinesY - point2);
                            switchlabels2 = (comparisonsY == nr_nullclinesY - point2);
                            if flip1==1
                                pairingsY(switchlabels) = point1  + nr_nullclinesY;
                                comparisonsY(switchlabels2) = point1  + nr_nullclinesY;
                            else
                                pairingsY(switchlabels) = point1;
                                comparisonsY(switchlabels2) = point1;
                            end
                        else 
                            % when it says end of 2 replace with end of 1
                            switchlabels = (pairingsY == point2 + nr_nullclinesY);
                            switchlabels2 = (comparisonsY == point2 + nr_nullclinesY);
                            if flip1==1
                                pairingsY(switchlabels) = point1 + nr_nullclinesY;
                                comparisonsY(switchlabels2) = point1 + nr_nullclinesY;
                            else
                                pairingsY(switchlabels) = point1;
                                comparisonsY(switchlabels2) = point1;
                            end

                        end
                    end                                                
                else
                    % remove this combo from pairs and comparisons
                    pairingsY(bestPairIdx(1),:) = [];
                    
                    row1 = [point1 point2];
                    row2 = [point2 point1];

                    idx1 = ismember(comparisonsY,row1,"row");
                    idx2 = ismember(comparisonsY,row2,"row");

                    comparisonsY(idx1 | idx2,:) = [];

                    anglesNear(bestPairIdx(1)) = [];

                end
            end

        end
        % delete extra nullclines
        YbranchesXdata(:,nullclinesToDeleteY) = [];
        YbranchesYdata(:,nullclinesToDeleteY) = [];

        nr_nullclinesY = size(YbranchesXdata,2);

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
        fimplicit(ax,NIO,[minx-PlotToleranceX maxx+PlotToleranceX miny-PlotToleranceY maxy+PlotToleranceY],'Color',colr,'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
        
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
        fimplicit(ax,NIO,[minx-PlotToleranceX maxx+PlotToleranceX miny-PlotToleranceY maxy+PlotToleranceY],'Color',colr,'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
        
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



    