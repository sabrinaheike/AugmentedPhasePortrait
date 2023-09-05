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

% Assumption: i) X-nullclines can be expressed in terms of functions in x;
% ii) Y-nullclines can be expressed in terms of functions in y;


% Note: Clear any x and y values you may have in your matlab
% memory by `clear x y' then type `syms x y' before calling
% the function

%%%%%%%%%%%%%%%%%%%%EXAMPLE HOW TO RUN THE CODE
% Type in the Matlab console:
% > clear x y
% > syms x y
% To get the augmented phase portrait for the system: 2*x/(1+x+0.3*y),
% 3*y/(1+2*y+0.6*x) plotted in [0,2]x[0,3] with default values, type:

% > augmented_xy(2*x/(1+x+0.3*y), 3*y/(1+2*y+0.6*x),0, 2,0, 3)
% [alternatively:]
% > clear x y
% > syms x y
% > f=2*x/(1+x+0.3*y);
% > g=3*y/(1+2*y+0.6*x);
% > augmented_xy(f,g,0,2,0,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[ ] = augmented_xy(fin,gin, minx, maxx, miny, maxy, varargin)
    syms x y
 
    trivh = 0;
    trivk = 0;
    counth = 0;
    countk = 0;
    gvax0 = 0;
    gvay0 = 0;
    fvax0 = 0;
    fvay0 = 0;
    
    
    %%%%%%%%%%%%%%%%1: identify trivial nullclines and embedd functions in functions of x and y
    [Nf, Df] = numden(fin);
    [Ng, Dg] = numden(gin);
    
    %%%%%%%%%%%%%1.1: For f(x,y):
    if length(symvar(fin))==1
        % if f is only a function of one of the two variables, either x
        % or y, then embedd it into a function of x and y
        if symvar(fin)==x
            disp('The X-equation does not contain an Y-value, so it cannot be solved for Y. Please use a different code.');
            return
        else 
            f1=fin+abs(x-minx)-(x-minx);
        end
    else
            f1 = fin;
    end
    
    
    if subs(Df,x,0)==0
        % check if f(0,y) is well-defined
        disp('the Y-axis is a vertical asymptote for the X-equation');
        fvax0 = 1;    
    elseif subs(fin,x,0)==0
        % case when f(0,y)=0-> identifies trivial X-nullcline as Y-axis
        % (X=0)
        disp('the Y-axis is a trivial nullcine for the X-equation');
        trivh = 1;
    end
    
    
    if subs(Df,y,0)==0
        % check if f(x,0) is well-defined
        disp('the X-axis is a vertical asymptote for the X-equation');
        fvay0 = 1;
    elseif subs(fin,y,0)==x
        % case when f(x,0)=x-> identifies trivial X-nullcline as X-axis
        % (Y=0)
        disp('the X-axis is a trivial nullcine for the X-equation');
        counth = 1;
    end

    

    %%%%%%%%%%%%%1.2: For g(x,y):
    if length(symvar(gin))==1 
        % if g is only a function of one of the two variables, either x
        % or y, then embedd it into a function of x and y
        if symvar(gin)==y
            disp('The Y-equation does not contain an X-value, so it cannot be solved for X. Please use a different code.');
            return
        else   
            g1 = gin+abs(y-miny)-(y-miny);
        end
    else 
        g1 = gin;
    end
        
    if subs(Dg,x,0)==0
        % check if g(0,y) is well-defined
        disp('the Y-axis is a vertical asymptote for the Y-equation');
        gvax0 = 1;
    elseif subs(gin,x,0)==y
        % case when g(0,y)=y-> identifies trivial Y-nullcline as Y-axis
        % (X=0)
        disp('the Y-axis is a trivial nullcine for the Y-equation');
        trivk = 1;
    end

    
    if subs(Dg,y,0)==0
        % check if g(x,0) is well-defined
        disp('the X-axis is a vertical asymptote for the Y-equation');
        gvay0 = 1;
    elseif subs(gin,y ,0)==0
        % case when g(x,0)=0-> identifies trivial Y-nullcline as X-axis
        % (Y=0)
        disp('the X-axis is a trivial nullcine for the Y-equation');
        countk = 1;
    end
    
    
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%2: Check if optional parameters were provided
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%3: Calculate the nullclines
    H = solve(fin-x==0,y, 'Real', true); %X-nullcline 
    disp('The X-nullclines are (as functions in X)');
    disp(H);
    K = solve(gin-y==0,x, 'Real', true); %Y-nullcline
    disp('The Y-nullclines are (as functions in Y)');
    disp(K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%3.2 Check which nullclines are nontrivial and in region
    x_toplot = [];
    if size(H,1)-counth>0
        x_idx = false(size(H,1),1);      
        for i=1:size(H,1)
            plotrh = 1; %boolean to check if nullcline is nontrivial and within diagram
            h = H(i);
            
            
            hd = diff(h,x);
            [Nhd, Dhd] = numden(hd);
            va = solve(Dhd==0,x,  'Real', true);
            vax = va(va>minx);
            vax = vax(vax<maxx);
            if length(vax)>0
                disp('take care there is a vertical asymptote of this nullcline but I do not distinguish colors');
                disp(h);
                disp('you might be better off trying to express h as a function of y');
            end 

            if isnumeric(eval(h))==1 && eval(h)==0 
                plotrh = 0; %trivial nullcline
                disp('Be reminded: I am not including signs of next-iterate operators for trivial nullclines');
            end
            
            %check if the nullcline y=h(x) is within the desired region
            if symvar(h)==x
                hf = matlabFunction(h);
                xtest = minx:0.01:maxx;
                countpos = 0;
                for tx = 1:length(xtest)
                    xt1 = xtest(tx);
                    ht = feval(hf,xt1);
                    if ht>miny && ht<maxy
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
        x_toplot = H(x_idx);
    end

    y_toplot = [];
    if size(K,1)-trivk>0
        y_idx = false(size(K,1),1); 
        for j=1:size(K,1)
            plotrk = 1;% boolean to check if nullcline is nontrivial and within diagram
            k = K(j);
            
            kd = diff(k,y);
            [Nkd, Dkd] = numden(kd);
            wa = solve(Dkd==0,y,  'Real', true);
            way = wa(wa>miny);
            way = way(way<maxy);
            if length(way)>0
                disp('take care there is a vertical asymptote of this nullcline but I do not distinguish colors');
                disp(k);
                disp('you might be better off trying to express k as a function of x');
            end 
         
            
            if isnumeric(eval(k))==1 && eval(k)==0 
                plotrk = 0; %trivial nullcline
                disp('Be reminded: I am not including signs of next-iterate operators for trivial nullclines');
            end
            
            %check if the nullcline x=k(y) is within the desired region
            if symvar(k)==y
                kf = matlabFunction(k);
                ytest = miny:0.01:maxy;
                countposg = 0;
                for ty = 1:length(ytest)
                    yt1 = ytest(ty);
                    kt = feval(kf,yt1);
                    if kt>minx && kt<maxx
                        countposg = 1;
                    end
                end
                
                if countposg == 0
                    plotrk = 0;
                end
            end
            
            if plotrk ==1
                 y_idx(j) = 1;
            end
        end
        y_toplot = K(y_idx);
     end

    %%%%%%%%%%%%%%%%4: Set up grid
    
    % create the coordinate-system with step-size for evaluations
    xaxis = minx:(maxx-minx)/acc:maxx;  
    yaxis = miny:(maxy-miny)/acc:maxy;
    nr = max(1,max(length(x_toplot),length(y_toplot)));
    xaxish = xaxis(1:nr+2:end);
    yaxish = yaxis(1:nr+2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    %%%%%%%%%%%%%%%%%%%%%%%5: DRAW THE DIRECTION FIELD
    drawArrow = @(w,z, varargin) quiver( w(1),z(1),w(2)-w(1),z(2)-z(1),2, 'MaxHeadSize', 12, varargin{:});    
    ff = matlabFunction(f1);
    gf = matlabFunction(g1);
    for md=1:length(xaxish)-1
        for nd = 1:length(yaxish)-1
       % points to be evaluated
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
   

    %%%%%%%%%%%%%%%%6: Plot X root-curves, and signs of
    %%%%%%%%%%%%%%%  corresponding next-iterate operators
   
    symbolTolerancex = 0.1/(maxx-minx);

    % Set colours
    if length(x_toplot) <= 3
        XColours = [0.3, 0.45, 1; 0, 0, 1; 0, 0, 0.6];
    else
        XColours = colormap(winter(length(x_toplot)+4));
    end    

        %for nontrivial nullclines within region of interest
     for i = 1:length(x_toplot)
         h = x_toplot(i);                
                
        %a) Construct the associated next-iterate operator 
        if  length(symvar(h))==0
            %if nullcline is constant
            Lh = gin-h;
        else
            Lh = gin-subs(h, x, fin);  
        end
               
        %b) Plot the associated root-curve
            %make sure it is a function of x and y
        if length(symvar(Lh))<2
            if symvar(Lh)==y
                Lh = Lh+abs(x-minx)-(x-minx);
            elseif symvar(Lh)==x
                Lh = Lh+abs(y-miny)-(y-miny);
            else 
                Lh = Lh+abs(y-miny)-(y-miny)+abs(x-minx)-(x-minx);
            end
        end
        hold on
        fimplicit(Lh,'Color',XColours(i,:),'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
               
               
        %c) Plot the sign of this next-iterate operator 
        xaxish = xaxis(i+1:nr+2:end);
        yaxish = yaxis(i+2:nr+2:end);
        Lhf = matlabFunction(Lh);
        if drawred==1
            for m=1:length(xaxish)
                for n = 1:length(yaxish)
                    xt = xaxish(m);
                    yt = yaxish(n);
                    Lhv = feval(Lhf,xt, yt);
                    if isreal(Lhv)==1 && Lhv >symbolTolerancex &&xt>minx &&yt>miny && ff(xt,yt)>varargin{2} && gf(xt,yt)>varargin{3}
                        hold on
                        plot([xt,xt],[yt,yt],'MarkerFaceColor',XColours(i,:),'MarkerEdgeColor',XColours(i,:),'MarkerSize',8,...
                                'Marker','+',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    elseif isreal(Lhv)==1 && Lhv< -symbolTolerancex &&xt>minx &&yt>miny && ff(xt,yt)>varargin{2} && gf(xt,yt)>varargin{3}
                        hold on
                        plot([xt,xt],[yt,yt],'MarkerFaceColor',XColours(i,:),'MarkerEdgeColor',XColours(i,:),'MarkerSize',10,...
                                'Marker','_',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    end      
                end
            end
        else
            for m=1:length(xaxish)
                for n = 1:length(yaxish)
                    xt = xaxish(m);
                    yt = yaxish(n);
                    Lhv = feval(Lhf,xt, yt);
                    if isreal(Lhv)==1 && Lhv >symbolTolerancex &&xt>minx &&yt>miny
                        hold on
                        plot([xt,xt],[yt,yt],'MarkerFaceColor',XColours(i,:),'MarkerEdgeColor',XColours(i,:),'MarkerSize',8,...
                                'Marker','+',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    elseif isreal(Lhv)==1 && Lhv< -symbolTolerancex &&xt>minx &&yt>miny
                        hold on
                        plot([xt,xt],[yt,yt],'MarkerFaceColor',XColours(i,:),'MarkerEdgeColor',XColours(i,:),'MarkerSize',10,...
                                'Marker','_',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    end      
                end
            end
        end
     end

  %%%%%%%%%%%%%%%%7: Plot Y root-curves, and signs of
  %%%%%%%%%%%%%%%  corresponding next-iterate operators
  
    symbolTolerancey = 0.1/(maxy-miny);

    % Set colours
    if length(y_toplot) <= 3
        YColours = [0.8, 0, 0; 1, 0, 0; 0.5, 0, 0];
    else
        YColours = colormap(autumn(length(y_toplot)+4));
    end    

        % if nontrivial and within the region of interest, then
    for j = 1:length(y_toplot)
        k = y_toplot(j);
     
        %a) Construct the associated next-iterate operator 
        if  length(symvar(k))==0
            %if constant nullcline
            Lk = fin-k;
        else
            Lk = fin-subs(k, y, gin);  
        end
               
        %b) Plot the associated root-curve 
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
        hold on
        fimplicit(Lk,'Color',YColours(j,:),'MeshDensity', 500,'LineWidth',1.5,'HandleVisibility','off');
               
        %c) Plot the sign of this next-iterate operator 
        xaxish = xaxis(j+2:nr+2:end);
        yaxish = yaxis(j+1:nr+2:end);
        Lkf = matlabFunction(Lk);
        if drawred==1
            for m=1:length(xaxish)
                for n = 1:length(yaxish)
                    xt = xaxish(m);
                    yt = yaxish(n);
                    Lkv = feval(Lkf,xt, yt);
                    if isreal(Lkv)==1 && Lkv >symbolTolerancey &&xt>minx &&yt>miny && ff(xt,yt)>varargin{2} && gf(xt,yt)>varargin{3}
                        hold on
                        plot([xt,xt],[yt,yt],'MarkerFaceColor',YColours(j,:),'MarkerEdgeColor',YColours(j,:),'MarkerSize',8,...
                                'Marker','+',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    elseif isreal(Lkv)==1 && Lkv< -symbolTolerancey &&xt>minx &&yt>miny && ff(xt,yt)>varargin{2} && gf(xt,yt)>varargin{3}
                    hold on
                    plot([xt,xt],[yt,yt],'MarkerFaceColor',YColours(j,:),'MarkerEdgeColor',YColours(j,:),'MarkerSize',10,...
                                'Marker','_',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    end      
                end
            end
        else
            for m=1:length(xaxish)
                for n = 1:length(yaxish)
                    xt = xaxish(m);
                    yt = yaxish(n);
                    Lkv = feval(Lkf,xt, yt);
                    if isreal(Lkv)==1 && Lkv >symbolTolerancey &&xt>minx &&yt>miny
                        hold on
                        plot([xt,xt],[yt,yt],'MarkerFaceColor',YColours(j,:),'MarkerEdgeColor',YColours(j,:),'MarkerSize',8,...
                                'Marker','+',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    elseif isreal(Lkv)==1 && Lkv< -symbolTolerancey &&xt>minx &&yt>miny
                    hold on
                    plot([xt,xt],[yt,yt],'MarkerFaceColor',YColours(j,:),'MarkerEdgeColor',YColours(j,:),'MarkerSize',10,...
                                'Marker','_',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                    end      
                end
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
        plot([minx maxx],[0 0],'--', 'Color', [0.5, 0.8, 1], 'LineWidth', 3, 'HandleVisibility','off');
    end
    if trivk>0
        hold on
        plot([0 0],[miny maxy], '--','Color', [1, 0.5, 0.8], 'LineWidth', 2.5, 'HandleVisibility','off');
    end
    if countk>0
        hold on
        plot([minx maxx],[0 0],'--', 'Color', [1, 0.25, 0.4], 'LineWidth', 2.5, 'HandleVisibility','off');
    end

    % Non-trivial X-nullclines
    for i = 1:length(x_toplot)
         h = x_toplot(i);
         if isnumeric(eval(h))==1
             % constant nullcline
             hold on
             plot([minx,maxx],[eval(h) eval(h)], 'LineStyle', '--', 'Color', XColours(i,:), 'LineWidth', 3,'DisplayName','X-Nullcline');
         else
             hf = matlabFunction(h);
             hold on
             fplot(hf, [minx,maxx],'LineStyle', '--','Color', XColours(i,:),'LineWidth', 3,'DisplayName','X-Nullcline'); %;,'LineStyle', '--','Color', CMH(i+counth,:),'LineWidth', 2,'DisplayName','X-Nullcline'); 
         end
    end
    
    % Non-trivial Y-nullclines
    for j = 1:length(y_toplot)
        k = y_toplot(j);
        if isnumeric(eval(k))==1
            %constant nullcline
            hold on
            plot([eval(k) eval(k)],[miny, maxy], 'LineStyle', '--', 'Color', YColours(j,:), 'LineWidth', 2.5,'DisplayName','Y-Nullcline');
        else
            kf = matlabFunction(k);
            hold on
            fimplicit(kf-x,'LineStyle', '--', 'Color', YColours(j,:),'MeshDensity', 500, 'LineWidth', 2.5,'DisplayName','Y-Nullcline');
        end
    end
     
    
    %%%%%%%%%%%9: ADD EQUILIBRIA IN THE PLOT
    %a) equilibria due to intersection of nontrivial nullclines
    for e1 = 1:size(H,1)
        h = H(e1);
        for e2 = 1:size(K,1)
            k = K(e2);
            try
                lastwarn('')
                kcomph = subs(k,y,h);
                xe = solve(x==kcomph,x,'Real',true);
                    for xes=1:size(xe,1)
                       xx = xe(xes);
                       xSubx = subs(Df,x,xx);
                       ySubx = subs(Dg,x,xx);
                       if xSubx ~= 0 && ySubx ~=0
                           ye = subs(h,x,eval(xx));
                           if subs(xSubx,y,ye) ~= 0 && subs(ySubx,y,ye) ~= 0
                               hold on
                               plot([eval(xx),eval(xx)],[eval(ye),eval(ye)],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                                'Marker','o','LineStyle','none', 'HandleVisibility','off');
                           end
                       end
                    end
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    h2 = h == y;
                    k2 = k == x;
                    equilibria(h2,k2,minx,maxx,miny,maxy,Df,Dg)
                end
            catch
                h2 = h == y;
                k2 = k == x;
                equilibria(h2,k2,minx,maxx,miny,maxy,Df,Dg)
            end
        end
    end
    
    %b) covers again the y-axis, i.e. x=0
    if trivh>0 && gvax0==0
        if trivk>0
            disp("All points on the y-axis are equilibria")
        else
            for e2 = 1:size(K,1)
                k = K(e2);
                lastwarn('')
                y0 = solve(k==0,y,'Real',true);
                for yes=1:size(y0,1)
                    ye0 = y0(yes);
                    if subs(Df,[x,y],[0,ye0]) ~= 0 && subs(Dg,[x,y],[0,ye0]) ~= 0 
                        hold on
                        plot([0,0],[eval(ye0),eval(ye0)],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                        'Marker','o','LineStyle','none', 'HandleVisibility','off');
                    end
                end
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    h_triv = x == 0;
                    k2 = k == x;
                    equilibria(h_triv,k2,minx,maxx,miny,maxy,Df,Dg)
                end
            end
        end
    end
    
    %c) equilibria at trivial nullclines [covers the x-axis , i.e y=0]
    if countk>0  && fvay0==0
        if counth>0
            disp("All points on the x-axis are equilibria")
        else
            for e3 = 1:size(H,1)
                h = H(e3);
                lastwarn('')
                x0 = solve(h==0,x,'Real',true);
                for xes=1:size(x0,1)
                    xe0 = x0(xes);
                    if subs(Df,[x,y],[xe0,0]) ~= 0 && subs(Dg,[x,y],[xe0,0]) ~= 0 
                        hold on
                        plot([xe0,xe0],[0,0],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                        'Marker','o','LineStyle','none', 'HandleVisibility','off');
                    end
                end
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    h2 = h == y;
                    k_triv = y == 0;
                    equilibria(h2,k_triv,minx,maxx,miny,maxy,Df,Dg) %method 2
                end
            end
        end
    end
      
%     
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