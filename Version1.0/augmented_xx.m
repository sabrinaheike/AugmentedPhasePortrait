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

% > augmented_xx(2*x/(1+x+0.3*y), 3*y/(1+2*y+0.6*x),0, 2,0, 3)
% [alternatively:]
% > clear x y
% > syms x y
% > f=2*x/(1+x+0.3*y);
% > g=3*y/(1+2*y+0.6*x);
% > augmented_xx(f,g,0,2,0,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[ ] = augmented_xx(fin,gin, minx, maxx, miny, maxy, varargin)
    syms x y
    CMHo = colormap(summer(20));
    CMKo = colormap(winter(20));
 
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
        hold on 
        plot([0 0],[miny maxy],'--', 'Color', CMHo(1,:), 'LineWidth', 2, 'HandleVisibility','off');
        disp('the Y-axis is a trivial nullcine for the X-equation');
        trivh = 1;
    end
    

    if subs(Df,y,0)==0
        % check if f(x,0) is well-defined
        disp('the X-axis is a vertical asymptote for the X-equation');
        fvay0 = 1;
    elseif subs(fin,y,0)==x
        % case when f(x,y)=x-> identifies trivial X-nullcline as X-axis
        % (Y=0)
        hold on 
        plot( [minx maxx],[0 0],'--', 'Color', CMHo(1,:), 'LineWidth', 2, 'HandleVisibility','off');
        disp('the X-axis is a trivial nullcine for the X-equation');
        counth = 1;
    end

    

    %%%%%%%%%%%%%1.2: For g(x,y):
    if length(symvar(gin))==1 
        % if g is only a function of one of the two variables, either x
        % or y, then embedd it into a function of x and y
        if symvar(gin)==y
            g1 = gin+abs(x-minx)-(x-minx);
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
        % case when g(0,y)=y -> identifies the Y-axis as nontrivial
        % Y-nullcline
        hold on 
        plot( [0 0],[miny maxy], '--','Color', CMKo(1,:), 'LineWidth', 2, 'HandleVisibility','off');
        disp('the Y-axis is a trivial nullcine for the Y-equation');
        trivk = 1;
    end

    
    if subs(Dg,y,0)==0
        % check if g(x,0) is well-defined
        disp('the X-axis is a vertical asymptote for the Y-equation');
        gvay0 = 1;
    elseif subs(gin,y ,0)==0
        % case when g(x,0)=0 -> identifies the X-axis as trivial nullcline
        % for Y-equation
        hold on 
        plot( [minx maxx],[0 0],'--', 'Color', CMKo(1,:), 'LineWidth', 2, 'HandleVisibility','off');
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
    K = solve(gin-y==0,y, 'Real', true); %Y-nullcline
    disp('The Y-nullclines are (as functions in X)');
    disp(K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
    %%%%%%%%%%%%%%%%4: Set up: Grid and colors for root-curves, and signs
    %%%%%%%%%%%%%%%%of next-iterate operators
    if size(H,1)+4<20
        CMH = CMHo;
    else
        CMH= colormap(summer(size(H,1)+4));
    end
    if size(K,1)+4<20
        CMK = CMKo;
    else    
        CMK= colormap(winter(size(K,1)+4));
    end
    
    % create the coordinate-system with step-size for evaluations
    xaxis = minx:(maxx-minx)/acc:maxx;  
    yaxis = miny:(maxy-miny)/acc:maxy;
    nr = max(1,size(H,1)+size(K,1)-counth - countk);
    xaxish = xaxis(1:nr+1:end);
    yaxish = yaxis(1:nr+1:end);
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
   

    %%%%%%%%%%%%%%%%6: Plot (nontrivial) X-nullcines (in the interval), root-curves, and signs of
    %%%%%%%%%%%%%%%  corresponding next-iterate operators
    if size(H,1)-counth>0
        for i=1:size(H,1)
            plotrh = 1; %boolean to check if nullcline is nontrivial and within diagram
            h = H(i);
            hf = matlabFunction(h); 
            
            % if vertical asymptote in nullcline tell that
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
                % if nontrivial nullcline and within region
                %a) plot the nullcline
                if isnumeric(eval(h))==1
                   %constant nullcline
                    hold on
                    plot([minx,maxx],[eval(h) eval(h)], 'LineStyle', '--', 'Color', CMH(i*2+counth,:), 'LineWidth', 2,'DisplayName','X-Nullcline');
                else
                    hold on
                    fplot(hf, [minx,maxx],'LineStyle', '--', 'Color', CMH(i*2+counth,:), 'LineWidth', 2,'DisplayName','X-Nullcline');
                end
     
               
                %b) Construct the associated next-iterate operator
                if  length(symvar(h))==0
                    % constant nullcline
                    Lh = gin-h;
                else
                    Lh = gin-subs(h, x, fin);  
                end
               
               
                %c) Plot the associated root-curve
                % make sure it is a function in x and y
                if length(symvar(Lh))<2
                   if symsvar(Lh)==y
                        Lh = Lh+abs(x-minx)-(x-minx);
                   elseif symvar(Lh)==x
                        Lh = Lh+abs(y-miny)-(y-miny);
                   else 
                        Lh = Lh+abs(y-miny)-(y-miny)+abs(x-minx)-(x-minx);
                   end
               end
               hold on
               fimplicit(Lh,'Color',CMH(i*2+counth,:),'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
               
               %d) Plot the sign of this next-iterate operator 
               
               xaxish = xaxis(i+1:nr+1:end);
               yaxish = yaxis(i+1:nr+1:end);
               Lhf = matlabFunction(Lh);
               for m=1:length(xaxish)
                   for n = 1:length(yaxish)
                        xt = xaxish(m);
                        yt = yaxish(n);
                        Lhv = feval(Lhf,xt, yt);
                        if isreal(Lhv)==1 && Lhv >0 &&xt>minx &&yt>miny
                            hold on
                            plot([xt,xt],[yt,yt],'MarkerFaceColor',CMH(i*2+counth,:),'MarkerEdgeColor',CMH(i*2+counth,:),'MarkerSize',8,...
                                'Marker','+',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                        elseif isreal(Lhv)==1 && Lhv<0 &&xt>minx &&yt>miny
                            hold on
                            plot([xt,xt],[yt,yt],'MarkerFaceColor',CMH(i*2+counth,:),'MarkerEdgeColor',CMH(i*2+counth,:),'MarkerSize',10,...
                                'Marker','_',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                        end      
                   end
               end
            end
        end
    end
  

    %%%%%%%%%%%%%%%%7: Plot (nontrivial) Y-nullcines (in the interval), root-curves, and signs of
    %%%%%%%%%%%%%%%  corresponding next-iterate operators
    

        
    if size(K,1)-countk>0
        
        for j=1:size(K,1)
            plotrk = 1; % boolean to check if nullcline is nontrivial and within diagram
            k = K(j);
            
            % if there is a vertical asymptote of the nullcline -> at least
            % tell the user
            kd = diff(k,x);
            [Nkd, Dkd] = numden(kd);
            wa = solve(Dkd==0,x,  'Real', true);
            wax = wa(wa>minx);
            wax = wax(wax<maxx);
            if length(wax)>0
                disp('take care there is a vertical asymptote of this nullcline but I do not distinguish colors');
                disp(k);
                disp('you might be better off trying to express k as a function of y');
            end 
            
            
            kf = matlabFunction(k); 
            if isnumeric(eval(k))==1 && eval(k)==0 
                plotrk = 0; %trivial nullcline
                disp('Be reminded: I am not including signs of next-iterate operators for trivial nullclines');
            end
            %check if the nullcline y=k(x) is within the desired region
            if symvar(k)==x
                kf = matlabFunction(k);
                xtest2 = minx:0.01:maxx;
                countposk = 0;
                for tx2 = 1:length(xtest2)
                    xt2 = xtest2(tx2);
                    kt2 = feval(kf,xt2);
                    if kt2>miny && kt2<maxy
                        countposk = 1;
                    end
                end
                if countposk == 0
                    plotrk = 0;
                end
            end
            
            if plotrk ==1
                % if nontrivial nullcline and intersects the region of
                % interest, then 
                % a) plot the nullcline
                if isnumeric(eval(k))==1
                    hold on
                    plot([minx,maxx],[eval(k) eval(k)],'LineStyle', '--', 'Color', CMK(j*2+countk,:), 'LineWidth', 2,'DisplayName','Y-Nullcline');
                else
                    hold on
                    fplot(kf, [minx,maxx],'LineStyle', '--', 'Color', CMK(j*2+countk,:), 'LineWidth', 2,'DisplayName','Y-Nullcline');
                end
     
               
                %b) create the associated next-iterate operator 
                if  length(symvar(k))==0
                    Lk = gin-k;
                else
                    Lk = gin-subs(k, x, fin);  
                end
               
                %c) Plot the associated root-curve; make sure it is a
                %function of x and y
                if length(symvar(Lk))<2
                    if symvar(Lk)==y
                        Lk = Lk+abs(x-minx)-(x-minx);
                    elseif symvar(Lk)==x
                        Lk = Lk+abs(y-miny)-(y-miny);
                    else 
                        Lk = Lk+abs(y-miny)-(y-miny)+abs(x-minx)-(x-minx);
                    end
                end
               % plot the associated root-curve
               hold on
               fimplicit(Lk,'Color',CMK(j*2+countk,:),'LineWidth',1.5,'MeshDensity', 500,'HandleVisibility','off');
               
               
               %d)  Plot the sign of this next-iterate operator                
               xaxish = xaxis(j+1+size(H,1)-counth:nr+1:end);
               yaxish = yaxis(j:nr*1.5:end);
               Lkf = matlabFunction(Lk);
               for m=1:length(xaxish)
                   for n = 1:length(yaxish)
                        xt = xaxish(m);
                        yt = yaxish(n);
                        Lkv = feval(Lkf,xt, yt);
                        if isreal(Lkv)==1&&Lkv >0 &&xt>minx &&yt>miny
                            hold on
                            plot([xt,xt],[yt,yt],'MarkerFaceColor',CMK(j*2+countk,:),'MarkerEdgeColor',CMK(j*2+countk,:),'MarkerSize',8,...
                                'Marker','+',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                        elseif isreal(Lkv)==1&& Lkv<0 &&xt>minx &&yt>miny
                            hold on
                            plot([xt,xt],[yt,yt],'MarkerFaceColor',CMK(j*2+countk,:),'MarkerEdgeColor',CMK(j*2+countk,:),'MarkerSize',10,...
                                'Marker','_',...
                                'LineStyle','-','LineWidth',2, 'HandleVisibility','off');
                        end      
                   end
               end
            end
        end
    end

    %%%%%%%%%%%8: ADD EQUILIBRIA IN THE PLOT
    %a) equilibria due to intersection of nontrivial nullclines 
    
    for e1 = 1:size(H,1)
        h = H(e1);
        for e2 = 1:size(K,1)
            k = K(e2);
            xe = solve(h==k,x);
                for xes=1:size(xe,1)
                   xx = xe(xes);
                   ye = subs(h,x,eval(xx));
                    hold on
                    plot([eval(xx),eval(xx)],[eval(ye),eval(ye)],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                    'Marker','o','LineStyle','none', 'HandleVisibility','off');
                end
           
        end
    end
    
    %b) equilibria at trivial nullclines [covers the y-axis , i.e x=0]
    if trivh>0 && gvax0==0 
        for e2 = 1:size(K,1)
            k = K(e2);
            [Nk, Dk] = numden(k);
            if subs(Dk,x,0)==0
                return
            else
            y0 = subs(k,x,0);
                for yes=1:size(y0,1)
                      ye0 = y0(yes);
                      hold on
                      plot([0,0],[eval(ye0),eval(ye0)],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                        'Marker','o','LineStyle','none', 'HandleVisibility','off');
                end
            end
        end
    end
    
    %c) equilibria at trivial nullclines [covers again the y-axis , i.e x=0]
    if trivk>0  && fvax0==0
        for e3 = 1:size(H,1)
            h = H(e3);
            [Nh, Dh] = numden(h);
            if subs(Dh,x,0)==0
                continue
            else
                y0 = subs(h, x, 0);
                for yes=1:size(y0,1)
                      ye0 = y0(yes);
                      hold on
                      plot([0,0],[ye0,ye0],'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',10,...
                        'Marker','o','LineStyle','none', 'HandleVisibility','off');
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
        fimplicit(fcutf,'-.', 'Color','red', 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('f=', num2str(varargin{3})));
        gcut = gin-varargin{3};
        gcutf = matlabFunction(gcut);
        hold on
        fimplicit(gcutf, '-.', 'Color',[1 0.4 0], 'LineWidth', 2,'MeshDensity', 500, 'DisplayName', strcat('g=', num2str(varargin{2})));
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