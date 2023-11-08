using SymPy
using Plots

function augmented_xx(fin, gin, minx, maxx, miny, maxy; acc = 20, drawred = false, x_cutoff = nothing, y_cutoff = nothing)
    @syms x y

    trivh = 0
    trivk = 0
    counth = 0
    countk = 0
    gvax0 = 0
    gvay0 = 0
    fvax0 = 0
    fvay0 = 0

    # 1: identify trivial nullclines and embed functions in functions of x and y
    Nf, Df = numden(fin)
    Ng, Dg = numden(gin)

    # 1.1: For f(x, y):
    if length(symvar(fin)) == 1
        if symvar(fin) == x
            println("The X-equation does not contain a Y-value, so it cannot be solved for Y. Please use a different code.")
            return
        else
            f1 = fin + abs(x - minx) - (x - minx)
        end
    else
        f1 = fin
    end

    if subs(Df, x, 0) == 0
        println("The Y-axis is a vertical asymptote for the X-equation")
        fvax0 = 1
    elseif subs(fin, x, 0) == 0
        println("The Y-axis is a trivial nullcline for the X-equation")
        trivh = 1
    end

    if subs(Df, y, 0) == 0
        println("The X-axis is a vertical asymptote for the X-equation")
        fvay0 = 1
    elseif subs(fin, y, 0) == x
        println("The X-axis is a trivial nullcline for the X-equation")
        counth = 1
    end

    # 1.2: For g(x, y):
    if length(symvar(gin)) == 1
        if symvar(gin) == y
            g1 = gin + abs(x - minx) - (x - minx)
        else
            g1 = gin + abs(y - miny) - (y - miny)
        end
    else
        g1 = gin
    end

    if subs(Dg, x, 0) == 0
        println("The Y-axis is a vertical asymptote for the Y-equation")
        gvax0 = 1
    elseif subs(gin, x, 0) == y
        println("The Y-axis is a trivial nullcline for the Y-equation")
        trivk = 1
    end

    if subs(Dg, y, 0) == 0
        println("The X-axis is a vertical asymptote for the Y-equation")
        gvay0 = 1
    elseif subs(gin, y, 0) == 0
        println("The X-axis is a trivial nullcline for the Y-equation")
        countk = 1
    end

    # 2: Check if optional parameters were provided
    if drawred && (x_cutoff === nothing || y_cutoff === nothing)
        println("If you input a cutoff for X, you also need to input a cutoff for Y")
        return
    end

    # 3: Calculate the nullclines
    H = solve(fin - x, y, Real=true)
    println("The X-nullclines are (as functions in X):")
    for h in H
        println(h)
    end

    K = solve(gin - y, y, Real=true)
    println("The Y-nullclines are (as functions in X):")
    for k in K
        println(k)
    end

    # 3.2 Check which nullclines are nontrivial and in the region
    XnullclineInfo = []
    if length(H) - counth > 0
        x_idx = [false for _ in 1:length(H)]
        for i in 1:length(H)
            plotrh = true  # boolean to check if nullcline is nontrivial and within the diagram
            h = H[i]

            if isnumeric(eval(h)) && eval(h) == 0
                plotrh = false  # trivial nullcline
                println("Be reminded: I am not including signs of next-iterate operators for trivial nullclines")
            end

            # Check if the nullcline y=h(x) is within the desired region
            if symvar(h) == x
                hf = lambdify(h, [x])
                xtest = range(minx, stop=maxx, length=acc)
                ytest = hf(xtest)
                for j in 1:length(xtest)
                    if !isreal(ytest[j])  # ensure the nullcline is real
                        plotrh = false
                    end
                end
            end
            if plotrh
                push!(XnullclineInfo, [h, "r"])
            end
        end
    end
    if length(K) - countk > 0
        y_idx = [false for _ in 1:length(K)]
        for i in 1:length(K)
            plotrk = true  # boolean to check if nullcline is nontrivial and within the diagram
            k = K[i]

            if isnumeric(eval(k)) && eval(k) == 0
                plotrk = false  # trivial nullcline
                println("Be reminded: I am not including signs of next-iterate operators for trivial nullclines")
            end

            # Check if the nullcline x=k(y) is within the desired region
            if symvar(k) == y
                kf = lambdify(k, [y])
                ytest = range(miny, stop=maxy, length=acc)
                xtest = kf(ytest)
                for j in 1:length(ytest)
                    if !isreal(xtest[j])  # ensure the nullcline is real
                        plotrk = false
                    end
                end
            end
            if plotrk
                push!(XnullclineInfo, [k, "b"])
            end
        end
    end

    # 4: Plot the direction fields
    xvals = range(minx, stop=maxx, length=acc)
    yvals = range(miny, stop=maxy, length=acc)
    vxvals = [fin(x, y) for x in xvals, y in yvals]
    vyvals = [gin(x, y) for x in xvals, y in yvals]

    if x_cutoff !== nothing && y_cutoff !== nothing
        plot(legend=false, ratio=1, xlims=(minx, x_cutoff), ylims=(miny, y_cutoff))
    else
        plot(legend=false, ratio=1, xlims=(minx, maxx), ylims=(miny, maxy))
    end

    if gvax0 || gvay0 || fvax0 || fvay0
        println("Not generating quiver: One or more nullclines in diagram")
        println("Not generating quiver: A nullcline is a vertical asymptote in diagram")
        println("Not generating quiver: A nullcline is on the diagram")
        println("Not generating quiver: A nullcline is a trivial nullcline in diagram")
    else
        plot!(xvals, yvals, quiver=([vxvals vyvals]), color=:auto)
    end

    for i in 1:length(XnullclineInfo)
        nullcline = XnullclineInfo[i][1]
        color = XnullclineInfo[i][2]

        if x_cutoff === nothing || y_cutoff === nothing
            if symvar(nullcline) == x
                plot!(xvals, [lambdify(nullcline, [x => xx]) for xx in xvals], linecolor=color, linestyle=:dot)
            else
                yvals = range(miny, stop=maxy, length=acc)
                plot!([lambdify(nullcline, [y => yy]) for yy in yvals], yvals, linecolor=color, linestyle=:dot)
            end
        else
            if symvar(nullcline) == x
                plot!(xvals, [lambdify(nullcline, [x => xx]) for xx in xvals], linecolor=color, linestyle=:dot, xlims=(minx, x_cutoff), ylims=(miny, y_cutoff))
            else
                yvals = range(miny, stop=maxy, length=acc)
                plot!([lambdify(nullcline, [y => yy]) for yy in yvals], yvals, linecolor=color, linestyle=:dot, xlims=(minx, x_cutoff), ylims=(miny, y_cutoff))
            end
        end
    end

    if x_cutoff === nothing || y_cutoff === nothing
        xvals = range(minx, stop=maxx, length=acc)
        yvals = range(miny, stop=maxy, length=acc)
    else
        xvals = range(minx, stop=x_cutoff, length=acc)
        yvals = range(miny, stop=y_cutoff, length=acc)
    end

    # 5: Generate the equilibrium points
    equilibrium_points = []
    for xval in xvals, yval in yvals
        if Nf(x => xval, y => yval) == 0 && Ng(x => xval, y => yval) == 0
            push!(equilibrium_points, (xval, yval))
        end
    end

    println("Equilibrium Points:")
    for eq in equilibrium_points
        println(eq)
    end

    if isempty(equilibrium_points)
        println("No equilibrium points found.")
    end

    scatter!([p[1] for p in equilibrium_points], [p[2] for p in equilibrium_points], color=:red)

    return
end

# Example usage:
fin = x * (1 - x) - y
gin = x - y
minx = 0
maxx = 1
miny = 0
maxy = 1

augmented_xx(fin, gin, minx, maxx, miny, maxy, acc=10, drawred=true)
