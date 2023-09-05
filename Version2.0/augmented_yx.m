% Augmented phase portrait 
% Input: (f,g,minx, maxx, miny, maxy, acc^*, cutoffx^*, cutoffy^*)
% ^* indicates optional parameters
% 
% 
% f = functionhandle; right-hand side of the X-equation f = f(x,y)
% g = functionhandle; right-hand side of the Y-equation g = g(x,y) 
% minx = number; minimum x-value to be considered for plotting
% maxx = number; maximum x-value to be considered for plotting
% miny = number; minimum y-value to be considered for plotting
% maxy = number; maximum y-value to be considered for plotting

% optional parameters:
% acc = number; optional parameter for number of arrows and signs of next-iterate operators; default value = 20
% cutoffx  = number; if number given then a red dash-dotted curve is added in plot to show points (x,y)
% such that f(x,y)=cutoffx
% cutoffy =  number; if number given then a orange dash-dotted curve is added in plot to show points (x,y)
% such that g(x,y)=cutoffy

% Output: Plot of augmented phase portrait (dashed lines correspond to
% nullclines, solid lines to root-curves)

% Assumption: i) X-nullclines can be expressed in terms of functions in y;
% ii) Y-nullclines can be expressed in terms of functions in x;


% Note: Clear any x and y values you may have in your matlab
% memory by `clear x y' then type `syms x y' before calling
% the function

%%%%%%%%%%%%%%%%%%%%EXAMPLE HOW TO RUN THE CODE
% Type in the Matlab console:
% > clear x y
% > syms x y
% To get the augmented phase portrait for the system: 2*x/(1+x+0.3*y),
% 3*y/(1+2*y+0.6*x) plotted in [0,2]x[0,3] with default values, type:

% > augmented_yx(2*x/(1+x+0.3*y), 3*y/(1+2*y+0.6*x),0, 2,0, 3)
% [alternatively:]
% > clear x y
% > syms x y
% > f=2*x/(1+x+0.3*y);
% > g=3*y/(1+2*y+0.6*x);
% > augmented_yx(f,g,0,2,0,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[ ] = augmented_yx(fin,gin, minx, maxx, miny, maxy, varargin)
    syms x y
 
    trivh = 0;
    trivk = 0;
    counth = 0;
    countk = 0;
    gvay0 = 0;
    fvax0 = 0;

    

    %%1: identify trivial nullclines and embedd functions in functions
    %%of x and y
   
    [Nf, Df] = numden(fin);
    [Ng, Dg] = numden(gin);
    
    %%%%%%%%%%% 1.1: For f(x,y):
    if length(symvar(fin))==1 
        % if f is only a function of one of the two variables, either x
        % or y, then embedd it into a function of x and y
        if symvar(fin)==x
            f1=fin+abs(y-miny)-(y-miny);
        else 
            f1=fin+abs(x-minx)-(x-minx);
        end
    else
        % else, if f is a function of x and y, then no need to embedd
        f1 = fin;
    end
    
    if subs(Df,x,0)==0
            % check if f(0,y) is well-defined
            disp('the Y-axis is a vertical asymptote for the X-equation');
            fvax0 = 1;
    elseif subs(fin,x,0)==0
        % case when f(0,y)=0-> identifies trivial X-nullcline as Y-axis (X=0)
        disp('the Y-axis is a trivial nullcine for the X-equation');
        trivh = 1;
    end
    
    
    if subs(Df,y,0)==0
        % check if f(x,0) is well-defined
        disp('the X-axis is a vertical asymptote for the X-equation');
    elseif subs(fin,y,0)==x
        % case when f(x,0)=x-> identifies trivial X-nullcline as X-axis (Y=0)        
        disp('the X-axis is a trivial nullcine for the X-equation');
        counth = 1;
    end
    
    
    %%%%%%%%%%%%%%1.2: For g(x,y):
    if length(symvar(gin))==1 
        % if g is only a function of one of the two variables, either x
        % or y, then embedd it into a function of x and y
        if symvar(gin)==y
            g1 = gin+abs(x-minx)-(x-minx);
        else   
            g1 = gin+abs(y-miny)-(y-miny);
        end
    else 
       %else, g is already a function of x and y
        g1 = gin;
    end
        
    if subs(Dg,x,0)==0
        % check if g(0,y) is well-defined 
        disp('the Y-axis is a vertical asymptote for the Y-equation'); 
    elseif subs(gin,x,0)==y
        % case when g(0,y)=y-> identifies trivial Y-nullcline as Y-axis (X=0)    
        disp('the Y-axis is a trivial nullcine for the Y-equation');
        trivk = 1;
    end

    
    if subs(Dg,y,0)==0
        % check if g(x,0) is well-defined 
        disp('the X-axis is a vertical asymptote for the Y-equation');
        gvay0 = 1;
    elseif subs(gin,y ,0)==0
        % case when g(x,0)=0-> identifies trivial Y-nullcline as X-axis (Y=0)    
        disp('the X-axis is a trivial nullcine for the Y-equation');
        countk = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%2: check if optional parameters were provided
    drawred = 0;
    if nargin<7
        acc = 20;
    elseif nargin<8
        acc = varargin{1};
    elseif nargin<9
        disp('if you input cutoff for X, you also need to input a cutoff for Y');
        return 
    else 
        acc= varargin{1};
        drawred = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%3: Calculate the nullclines
    H = solve(fin-x==0,x, 'Real', true); %X-nullcline 
    disp('The X-nullclines are (as functions in Y)');
    disp(H);
    K = solve(gin-y==0,y, 'Real', true); %Y-nullcline
    disp('The Y-nullclines are (as functions in X)');
    disp(K);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%3.2 Check which nullclines are nontrivial and in region
    XnullclineInfo = {};
    if size(H,1)-trivh>0
        x_idx = false(size(H,1),1);
        for i=1:size(H,1)
            plotrh = 1; % boolean to check if nullcline is nontrivial and within diagram
            h = H(i);
                   
            if isnumeric(eval(h))==1 && eval(h)==0 
                plotrh = 0;  % trivial nullclines 
                disp('Be reminded: I am not including signs of next-iterate operators for trivial nullclines');
            end
            
            %check if the nullcline x=h(y) is within the desired region
            if symvar(h)==y
                hf = matlabFunction(h);
                ytest1 = miny:0.01:maxy;
                countpos = 0;
                for ty1 = 1:length(ytest1)
                    yt1 = ytest1(ty1);
                    ht = feval(hf,yt1);
                    if ht>minx && ht<maxx
                        countpos = 1;
                    end
                end
                if countpos == 0
                    plotrh = 0;
                end
            end
            
            if plotrh ==1
                x_idx(i) = 1;
            end
        end
        x_toplot = unique(H(x_idx)); % extract nullclines of interest discarding repeats
        XnullclineInfo = num2cell(x_toplot); %save list to cell array to collect info about nullclines
    end

    YnullclineInfo = {};
    if size(K,1)-countk>0
        y_idx = false(size(K,1),1);
        for j=1:size(K,1)
            plotrk = 1; % boolean to check if nullcline is nontrivial and within diagram
            k = K(j);
            
            if isnumeric(eval(k))==1 && eval(k)==0 
                plotrk = 0; %trivial nullcine
                disp('Be reminded: I am not including signs of next-iterate operators for trivial nullclines');
            end
            %check if the nullcline x=k(y) is within the desired region
            if symvar(k)==x
                kf = matlabFunction(k);
                xtest = minx:0.01:maxx;
                countposk = 0;
                for tx = 1:length(xtest)
                    xt = xtest(tx);
                    kt = feval(kf,xt);
                    if kt>miny && kt<maxy
                        countposk = 1;
                    end
                end
                if countposk == 0
                    plotrk = 0;
                end
            end
            
            if plotrk ==1
                y_idx(j) = 1;
            end
        end
        y_toplot = unique(K(y_idx));
        YnullclineInfo = num2cell(y_toplot);
    end
    
    %%%%%%%%%%%%%3.3 Check for asymptotes in nullclines
    global asymptoteList
    for i = 1:size(XnullclineInfo,1)
        h = XnullclineInfo{i,1};
        hd = diff(h,y);
        [~, Dhd] = numden(hd);
        lastwarn('')
        va = solve(Dhd==0,y,  'Real', true);
        asymptoteList = va(va>miny);
        asymptoteList = asymptoteList(asymptoteList<maxy);
        [warnMsg, ~] = lastwarn;
        if ~isempty(warnMsg) %if there was a warning using solve
            tolerance = 0.001;
            % try alternate method to solve for asymptotes
            asymptote(Dhd,miny+tolerance,maxy-tolerance,y)
        end
        if ~isempty(asymptoteList)
            disp('take care there is a vertical asymptote (grey dashed line) of this nullcline but I do not distinguish colors');
            disp(h);
            disp('you might be better off trying to express h as a function of x');
        end
        XnullclineInfo{i,2} = unique(asymptoteList);
    end
      
    for j = 1:size(YnullclineInfo,1)
        k = YnullclineInfo{j,1};    
        kd = diff(k,x);
        [~, Dkd] = numden(kd);
        lastwarn('')
        wa = solve(Dkd==0,x,  'Real', true);
        asymptoteList = wa(wa>minx);
        asymptoteList = asymptoteList(asymptoteList<maxx);
        [warnMsg, ~] = lastwarn;
        if ~isempty(warnMsg) %if there was a warning using solve
            tolerance = 0.001;
            % try alternate method to solve for asymptotes
            asymptote(Dkd,minx+tolerance,maxx-tolerance,x)
        end
        if ~isempty(asymptoteList)
            disp('take care there is a vertical asymptote (grey dashed line) of this nullcline but I do not distinguish colors');
            disp(k);
            disp('you might be better off trying to express k as a function of y');
        end
        YnullclineInfo{j,2} = unique(asymptoteList);
    end
  

    %%%%%%%%%%%%%%%%%4: Set up grid 
    
    % create the coordinate-system with step-size for evaluations
    xaxis = minx:(maxx-minx)/acc:maxx;  
    yaxis = miny:(maxy-miny)/acc:maxy;
    nr = max(1,max(size(XnullclineInfo,1),size(YnullclineInfo,1)));
    xaxish = xaxis(1:nr+2:end);  
    yaxish = yaxis(1:nr+2:end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    
    
    
    %%%%%%%%%%%%%%5:  DRAW THE DIRECTION FIELD
    drawArrow = @(w,z, varargin) quiver( w(1),z(1),w(2)-w(1),z(2)-z(1),2, 'MaxHeadSize', 12, varargin{:});    
    ff = matlabFunction(f1);
    gf = matlabFunction(g1);
    
    for md=1:length(xaxish)-1
        for nd = 1:length(yaxish)-1
            % point to be evaluated
            xt = xaxish(md);
            yt = yaxish(nd);
            fval = feval(ff, xt, yt);
            gval = feval(gf, xt, yt);
            
            % end of the arrowlength to be evaluated (to avoid arrows
            % pointing across nullclines in the wrong direction)
            arrowlex = (maxx-minx)/(2*acc); 
            arrowley = (maxy-miny)/(2*acc);
            xt2 = xt+sign(fval-xt)*2*arrowlex;
            yt2 = yt+sign(gval-yt)*2*arrowley; 
            fval2 = feval(ff, xt2, yt);
            gval2 = feval(gf, xt, yt2);
            % arrowlength
            if fval>xt && fval2>xt2 && gval>yt  && gval2>yt2
                hold on
                drawArrow([xt, xt+arrowlex],[yt, yt],'linewidth',1,'color','k','HandleVisibility','off');     
                hold on
                drawArrow([xt, xt],[yt, yt+arrowley], 'linewidth',1,'color','k','HandleVisibility','off');
            elseif fval>xt && fval2>xt2 && gval<yt  && gval2<yt2
                hold on
                drawArrow([xt, xt],[yt, yt-arrowley],'linewidth',1,'color','k','HandleVisibility','off');     
                hold on
                drawArrow([xt, xt+arrowlex],[yt, yt], 'linewidth',1,'color','k','HandleVisibility','off');

            elseif fval<xt  && fval2< xt2 && gval<yt  && gval2<yt2
                hold on
                drawArrow([xt, xt],[yt, yt-arrowley],'linewidth',1,'color','k','HandleVisibility','off');     
                hold on
                drawArrow([xt, xt-arrowlex],[yt, yt],'linewidth',1,'color','k','HandleVisibility','off');     
            elseif fval<xt  && fval2<xt2 && gval>yt  && gval2>yt2
                hold on
                drawArrow([xt, xt],[yt, yt+arrowley],'linewidth',1,'color','k','HandleVisibility','off');     
                hold on
                drawArrow([xt, xt-arrowlex],[yt, yt],'linewidth',1,'color','k','HandleVisibility','off');     
            end
        end
    end
   

      
    %%%%%%%%%%%%%%%%6.1: Construct Next Iterative Operator and Plot Root
    %%%%%%%%%%%%%%%%Curves for X-Equation

    % a) set colours
        
    % Keep track of how many colours needed due to different branches of
    % nullclines
    ColrAdjustX = 0;
    for i = 1:size(XnullclineInfo,1)
        XnullclineInfo{i,3} = i+ColrAdjustX:i+ColrAdjustX+length(XnullclineInfo{i,2});
        ColrAdjustX = ColrAdjustX + length(XnullclineInfo{i,2});
    end

    % Set colourscheme depending on total number of branches
    if ~isempty(XnullclineInfo)
        if XnullclineInfo{end,3}(end)<= 3
            XColours = [0.3, 0.45, 1; 0, 0, 1; 0, 0, 0.6];
        else
            XColours = colormap(winter(size(XnullclineInfo,1)+4));
        end
    end

    % Asign colours to each branch
    for i = 1:size(XnullclineInfo,1)
        XnullclineInfo{i,4} = XColours(XnullclineInfo{i,3},:);
    end   


        % if nontrivial and within region of interest, then 
    for i = 1:size(XnullclineInfo,1)
        h = XnullclineInfo{i,1};

        %b) Construct the associated next-iterate operator 
        if  length(symvar(h))==0
            % if nullcline is constant
            Lh = fin-h;
        else
            Lh = fin-subs(h, y, gin);  
        end


            % make sure it is a function of x and y
        if length(symvar(Lh))<2
            if symvar(Lh)==y
                Lh = Lh+abs(x-minx)-(x-minx);
            elseif symvar(Lh)==x
                Lh = Lh+abs(y-miny)-(y-miny);
            else 
                Lh = Lh+abs(y-miny)-(y-miny)+abs(x-minx)-(x-minx);
            end
        end

        Lhf = matlabFunction(Lh);
        XnullclineInfo{i,5} = Lhf;

        % c) plot root curve in average shade of nullcline branches
        hold on
        fimplicit(Lh,'Color',mean(XnullclineInfo{i,4},1),'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
    end


    %%%%%%%%%%%%%%%%6.2: Construct Next Iterative Operator and Plot Root
    %%%%%%%%%%%%%%%%Curves for Y-Equation

    % a) Set Colours    
        
    % Keep track of how many colours needed due to different branches of
    % nullclines
    ColrAdjustY = 0;
    for j = 1:size(YnullclineInfo,1)
        YnullclineInfo{j,3} = j+ColrAdjustY:j+ColrAdjustY+length(YnullclineInfo{j,2});
        ColrAdjustY = ColrAdjustY + length(YnullclineInfo{j,2});
    end

    % Set colourscheme depending on total number of branches  
    if ~isempty(YnullclineInfo)
        if YnullclineInfo{end,3}(end) <= 3
            YColours = [0.8, 0, 0; 1, 0, 0; 0.5, 0, 0];
        else
            YColours = colormap(autumn(size(YnullclineInfo,1)+4));
        end
    end

    % Asign colours to each branch
    for j = 1:size(YnullclineInfo,1)
        YnullclineInfo{j,4} = YColours(YnullclineInfo{j,3},:);
    end


        % if nontrivial and within the region of interest, then        
    for j = 1:size(YnullclineInfo,1)
        k = YnullclineInfo{j,1};
                
        %b) Construct the associated next-iterate operator 
        if  length(symvar(k))==0
            % if nullcline is constant   
            Lk = gin-k;
        else
            Lk = gin-subs(k, x, fin);  
        end
                
 
            % make sure it is a function of x and y
        if length(symvar(Lk))<2
            if symvar(Lk)==y
                Lk = Lk+abs(x-minx)-(x-minx);
            elseif symvar(Lk)==x
                Lk = Lk+abs(y-miny)-(y-miny);
            else 
                Lk = Lk+abs(y-miny)-(y-miny)+abs(x-minx)-(x-minx);
            end
        end

        Lkf = matlabFunction(Lk);
        YnullclineInfo{j,5} = Lkf;

        % c) plot the associated root-curve
        hold on
        fimplicit(Lk,'Color',mean(YnullclineInfo{j,4},1),'LineWidth',1.5, 'MeshDensity', 500,'HandleVisibility','off');
    end


    %%%%%%%%%%%%%%%%7: Plot signs of next-iterate operators

    PlotGrid = cell(acc+1,acc+1,3); % create grid to store symbol information
        % PlotGrid{m,n,1} = 1 if that point should be plotted
        % PlotGrid{m,n,2} contains the symbol to be plotted
        % PlotGrid{m,n,3} contains the colour of symbol

    %%%Evaluate sign of next iterative operator
        % X-Nullclines
    symbolToleranceX = 0.1/(maxx-minx); 
    for i = 1:size(XnullclineInfo,1) % for each nullcline
        for n=i+1:nr+2:length(xaxis)
            for m = i+2:nr+2:length(yaxis)
                xt = xaxis(n);
                yt = yaxis(m);
                Lhv = feval(XnullclineInfo{i,5},xt, yt); % evaluate NIO at point
                if isreal(Lhv)
                    if Lhv > symbolToleranceX % positive NIO
                        PlotGrid{m,n,1} = 1; % mark this point to be plotted
                        PlotGrid{m,n,2} = '+'; 
                    elseif Lhv < -symbolToleranceX % negative NIO
                        PlotGrid{m,n,1} = 1; % mark this point to be plotted
                        PlotGrid{m,n,2} = '_';
                    end
                end
            end
        end
    end

        % Y-Nullclines
    symbolToleranceY = 0.1/(maxy-miny);
    for j = 1:size(YnullclineInfo,1)
        for n=j+2:nr+2:length(xaxis)
            for m = j+1:nr+2:length(yaxis)
                xt = xaxis(n);
                yt = yaxis(m);
                Lkv = feval(YnullclineInfo{j,5},xt, yt);
                if isreal(Lkv)
                    if Lkv > symbolToleranceY
                        PlotGrid{m,n,1} = 1; % mark this point to be plotted
                        PlotGrid{m,n,2} = '+';
                    elseif Lkv < -symbolToleranceY
                        PlotGrid{m,n,1} = 1; % mark this point to be plotted
                        PlotGrid{m,n,2} = '_';
                    end
                end
            end
        end
    end
    
    %%% Make sure points are mapped above the cutoff values
    if drawred == 1
        for m=1:acc+1
            for n = 1:acc+1
                xt = xaxis(n);
                yt = yaxis(m);
                % if point is mapped below cutoff value for x or y 
                if ff(xt,yt)<varargin{2} || gf(xt,yt)<varargin{3}
                    PlotGrid{m,n,1} = 0; % mark it to NOT be plotted
                end
            end
        end
    end

    %%% Set colours
        % X-nullclines
    for i = 1:size(XnullclineInfo,1)
        if ~isempty(XnullclineInfo{i,2}) %if asymptote in nullcline
            for n=i+1:nr+2:length(xaxis)
                for m = i+2:nr+2:length(yaxis)
                    xt = xaxis(n);
                    yt = yaxis(m);
                    if gf(xt,yt) > XnullclineInfo{i,2}(end) %point mapped above all asymptotes
                        % set symbol colour to the colour of the last
                        % nullcline branch
                        PlotGrid{m,n,3} = XnullclineInfo{i,4}(end,:);
                    else
                        % if point gets mapped below r^th
                        % asymptote set colour to the same as the r^th
                        % nullcline branch
                        for r = 1:length(XnullclineInfo{i,2})
                            if gf(xt,yt) < XnullclineInfo{i,2}(r)
                                PlotGrid{m,n,3} = XnullclineInfo{i,4}(r,:);
                                break
                            end
                        end
                    end
                end
            end
        else
            % if no asymptote in nullcline set all colours of symbols to be
            % the same as the nullcline colour
            for n=i+1:nr+2:length(xaxis)
                for m = i+2:nr+2:length(yaxis)
                    PlotGrid{m,n,3} = XnullclineInfo{i,4};
                end
            end
        end
    end

        % Y-nullclines
    for j = 1:size(YnullclineInfo,1)
        if ~isempty(YnullclineInfo{j,2}) %if asymptote in nullcline
            for n=j+2:nr+2:length(xaxis)
                for m = j+1:nr+2:length(yaxis)
                    xt = xaxis(n);
                    yt = yaxis(m);
                    if ff(xt,yt) > YnullclineInfo{j,2}(end) %point mapped to the right of all asymptotes
                        % set symbol colour to the colour of the last
                        % nullcline branch
                        PlotGrid{m,n,3} = YnullclineInfo{j,4}(end,:);
                    else
                        % if point gets mapped to the left of r^th
                        % asymptote set colour to the same as the r^th
                        % nullcline branch
                        for r = 1:length(YnullclineInfo{j,2})                        
                            if ff(xt,yt) < YnullclineInfo{j,2}(r)
                                PlotGrid{m,n,3} = YnullclineInfo{j,4}(r,:);
                                break
                            end
                        end
                    end
                end
            end
        else 
            for n=j+2:nr+2:length(xaxis)
                for m = j+1:nr+2:length(yaxis)
                    PlotGrid{m,n,3} = YnullclineInfo{j,4};
                end
            end
        end
    end

    % plot all symbols
    for m = 1:acc+1
        for n = 1:acc+1
            if PlotGrid{m,n,1} == 1
                xt = xaxis(n);
                yt = yaxis(m);
                plot([xt,xt],[yt,yt],'MarkerFaceColor',PlotGrid{m,n,3},'MarkerEdgeColor',PlotGrid{m,n,3},'MarkerSize',8,...
                                'Marker',PlotGrid{m,n,2},...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
            end
        end
    end

    %%%%%%%%%%%8: Plot nullcines (in the interval)
    
    % Trivial nullclines
    if trivh>0
        hold on
        plot([0 0],[miny maxy],'--', 'Color', [0.5, 0.35,1], 'LineWidth', 3, 'HandleVisibility','off');
    end
    if counth>0
        hold on
        plot( [minx maxx],[0 0],'--', 'Color', [0.5, 0.8, 1], 'LineWidth', 3, 'HandleVisibility','off');
    end
    if trivk>0
        hold on
        plot( [0 0],[miny maxy], '--','Color', [1, 0.5, 0.8], 'LineWidth', 2.5, 'HandleVisibility','off');
    end
    if countk>0
        hold on
        plot( [minx maxx],[0 0],'--', 'Color', [1, 0.25, 0.4], 'LineWidth', 2.5, 'HandleVisibility','off');
    end

    % Non-trivial X-nullclines
    for i = 1:size(XnullclineInfo,1)
        h = XnullclineInfo{i,1};
        if isnumeric(eval(h))==1                   
            % constant nullcline
            hold on
            plot([eval(h) eval(h)],[miny,maxy], 'LineStyle', '--', 'Color', XnullclineInfo{i,4}, 'LineWidth', 3,'DisplayName','X-Nullcline');
        else
            hf = matlabFunction(h);
            yIntervals = [miny; XnullclineInfo{i,2}; maxy];
            plotTolerance = (maxy-miny)/10000;
            for m = 1:length(yIntervals) - 1 %plot each branch of nullcline in its assigned colour
                if yIntervals(m+1) - yIntervals(m) > 2*plotTolerance
                    hold on
                    fimplicit(hf-x, [minx maxx eval(yIntervals(m)+plotTolerance) eval(yIntervals(m+1)-plotTolerance)],'LineStyle', '--', 'Color', XnullclineInfo{i,4}(m,:), 'LineWidth', 3,'MeshDensity', 500,'DisplayName','X-Nullcline');
                end
            end
        end
    end

    % Non-Trivial Y-nullclines
    for j = 1:size(YnullclineInfo,1)
        k = YnullclineInfo{j,1};
        if isnumeric(eval(k))==1
            % constant nullcline
            hold on
            plot([minx, maxx],[eval(k) eval(k)], 'LineStyle', '--', 'Color', YnullclineInfo{j,4}, 'LineWidth', 2.5,'DisplayName','Y-Nullcline');
        else
            kf = matlabFunction(k);
            xIntervals = [minx; YnullclineInfo{j,2}; maxx];
            plotTolerance = (maxx-minx)/10000;
            for m = 1:length(xIntervals) - 1
                if xIntervals(m+1) - xIntervals(m) > 2*plotTolerance
                    hold on
                    fplot(kf, [eval(xIntervals(m)+plotTolerance),eval(xIntervals(m+1)-plotTolerance)],'LineStyle', '--', 'Color', YnullclineInfo{j,4}(m,:), 'LineWidth', 2.5,'DisplayName','Y-Nullcline');
                end
            end
        end
    end
    
    %%%%%%%%%%%9: Plot nullcline asymptotes
        % a) X-nullclines 
    for i=1:size(XnullclineInfo,1)
        for r = 1:length(XnullclineInfo{i,2})
            plot([minx maxx],[eval(XnullclineInfo{i,2}(r)), eval(XnullclineInfo{i,2}(r))], 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3, 'HandleVisibility','off');
        end
    end

        % b) Y-nullclines
    for j=1:size(YnullclineInfo,1)
        for r = 1:length(YnullclineInfo{j,2})
            plot([eval(YnullclineInfo{j,2}(r)), eval(YnullclineInfo{j,2}(r))], [miny maxy], 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3, 'HandleVisibility','off');
        end
    end

    %%%%%%%%%%%10: ADD EQUILIBRIA IN THE PLOT
    %a) equilibria due to intersection of nontrivial nullclines
    for e1 = 1:size(H,1)
        h = H(e1);
        for e2 = 1:size(K,1)
            k = K(e2);
            try
                lastwarn('')
                kcomph = subs(k,x,h);
                ye = solve(y==kcomph,y,'Real',true);
                for yes=1:size(ye,1)
                    yy = ye(yes);
                    xSuby = subs(Df,y,yy);
                    ySuby = subs(Dg,y,yy);
                    if xSuby ~= 0 && ySuby ~= 0
                        xe = subs(h,y,eval(yy));
                        if subs(xSuby,x,xe) ~= 0 && subs(ySuby,x,xe) ~= 0
                            hold on
                            plot([eval(xe),eval(xe)],[eval(yy),eval(yy)],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                            'Marker','o','LineStyle','none', 'HandleVisibility','off');
                        end
                    end
                end
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                     h2 = h == x;
                     k2 = k == y;
                     equilibria(h2,k2,minx,maxx,miny,maxy,Df,Dg) %method 2
                end
            catch
                h2 = h == x;
                k2 = k == y;
                equilibria(h2,k2,minx,maxx,miny,maxy,Df,Dg)  
            end
        end
    end
    
    %b) equilibria at trivial nullclines [covers the x-axis , i.e y=0]
    if counth>0 && gvay0==0
        if countk > 0
            disp("All points on the x-axis are equilibria")
        else
            for e2 = 1:size(K,1)
                k = K(e2);
                lastwarn('')
                x0 = solve(k==0,x,'Real',true);
                for xes=1:size(x0,1)
                    xe0 = x0(xes);
                    if subs(Df,[x,y],[xe0,0]) ~= 0 && subs(Dg,[x,y],[xe0,0]) ~= 0
                        hold on
                        plot([eval(xe0),eval(xe0)],[0,0], 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                        'Marker','o','LineStyle','none', 'HandleVisibility','off');
                    end
                end
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    h_triv = y == 0;
                    k2 = k == y;
                    equilibria(h_triv,k2,minx,maxx,miny,maxy,Df,Dg)
                end
            end
        end
    end
    
    %c) equilibria at trivial nullclines [covers the y-axis, i.e., x=0]
    if trivk>0  && fvax0==0
        if trivh>0
            disp("All points on the y-axis are equilibria")
        else
            for e3 = 1:size(H,1)
                h = H(e3);
                lastwarn('')
                y0 = solve(h==0,y,'Real',true);
                for yes=1:size(y0,1)
                    ye0 = y0(yes);
                    if subs(Df,[x,y],[0,ye0]) ~= 0 && subs(Dg,[x,y],[0,ye0]) ~= 0 
                        hold on
                        plot([0,0],[ye0,ye0],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                        'Marker','o','LineStyle','none', 'HandleVisibility','off');
                    end
                end
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    h2 = h == x;
                    k_triv = x == 0;
                    equilibria(h2,k_triv,minx,maxx,miny,maxy,Df,Dg) %method 2
                end
            end
        end
    end

    % intersection between 2 trivial nullclines
    if counth>0 && trivk>0
        if countk==0 && trivh==0
            if subs(Df,x,0) ~= 0 && subs(Dg,y,0) ~= 0
                plot([0,0],[0,0],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                            'Marker','o','LineStyle','none', 'HandleVisibility','off');
            end
        end
    end
    
    %%%%%%%%%%%%%%%%9:  Plot cutoff functions (if given)
    if drawred ==1
        fcut = fin-varargin{2};
        fcutf = matlabFunction(fcut);
        hold on
        fimplicit(fcutf,':', 'Color',[0.45 0.85 0.45], 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('f=', num2str(varargin{2})));
        gcut = gin-varargin{3};
        gcutf = matlabFunction(gcut);
        hold on
        fimplicit(gcutf, ':', 'Color',[0 0.45 0], 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('g=', num2str(varargin{3})));
    end
    
    
    %%%%%%%%%%%%%%%%%10: Set labels, size of diagram, and fontsize
    xlim([minx, maxx]);
    ylim([miny, maxy]);
    xlabel('X');
    ylabel('Y', 'rotation', 0);
    legend
    box on
    set(gca,'FontSize',18)
    hold off
end


% equilibria function to numerically find nullcline intersections when unable to solve symbolically 
function[] = equilibria(h_in,k_in,min_evalx,max_evalx,min_evaly,max_evaly,Denomf,Denomg)
    syms x y
    tolerance = 0.001;
    if max_evalx - min_evalx < tolerance || max_evaly - min_evaly < tolerance
        return
    end
    % stops after finding one intersection point 
    [x_eq,y_eq] = vpasolve(h_in,k_in, [x,y], [min_evalx, max_evalx; min_evaly, max_evaly]);
    % if one point was found, there may be more, so we need to exclude the
    % found point and search the remaining area by dividing it into regions
    % without the previously found point.
    if ~isempty(x_eq)
        if subs(Denomf,[x,y],[eval(x_eq(1)),eval(y_eq(1))]) ~=0 && subs(Denomg,[x,y],[eval(x_eq(1)),(y_eq(1))]) ~=0
            hold on
            plot(x_eq(1),y_eq(1),'k.','MarkerSize',40,'HandleVisibility','off')
        end
        if x_eq(1) - min_evalx > tolerance 
            % search area to the left of point
            equilibria(h_in,k_in,min_evalx,x_eq(1)-tolerance,min_evaly,max_evaly,Denomf,Denomg)
        end
        if max_evalx - x_eq(1) > tolerance 
            % search area to the right of point
            equilibria(h_in,k_in,x_eq(1)+tolerance,max_evalx,min_evaly,max_evaly,Denomf,Denomg)
        end
        if y_eq(1) - min_evaly > tolerance 
            % search area below point
            equilibria(h_in,k_in,x_eq(1)-tolerance,x_eq(1)+tolerance,min_evaly,y_eq(1)-tolerance,Denomf,Denomg)
        end
        if max_evaly - y_eq(1) > tolerance 
            % search area above point
            equilibria(h_in,k_in,x_eq(1)-tolerance,x_eq(1)+tolerance,y_eq(1)+tolerance,max_evaly,Denomf,Denomg) 
        end
    end
    hold off
end

% asymptote identification function to numerically find asymptotes when unable to solve symbolically 
function [] = asymptote(Denom,min_val,max_val,variable)
    %%% Alternate method to search for asymptotes in nullclines if solve
    %%% method fails
    tolerance = 0.001;
    global asymptoteList
    if max_val - min_val < tolerance
        return
    end
    % vpasolve returns only one solution unless the equation is polynomial
    xva = vpasolve(Denom==0,variable,[min_val, max_val]);
    % if one value was found, there may be more asymptotes, so we need to exclude the
    % found value and search the remaining range by dividing it into regions
    % without the previously found value.
    if ~isempty(xva)
        asymptoteList = [asymptoteList;xva(1)];
        if xva(1) - min_val > tolerance
            % search the area less than the found asymptote
            asymptote(Denom,min_val,xva(1)-tolerance,variable)
        end
        if max_val - xva(1) > tolerance
            % search the area greater than the found asymptote
            asymptote(Denom,xva(1)+tolerance,max_val,variable)
        end
    end
end
