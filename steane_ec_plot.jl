using Plots,CSV,DataFrames,ArgParse

function plot_mi_w_errbar()
    l=args["size"]
    r=args["rate"]

    plot()
    filename = "data/mi_L_$(l)_r_$(r).csv"
    df = CSV.read(filename, DataFrame)
    p_list = sort(unique(df[!,"p"]))
    # p_str_list = [replace(string(p),"."=>"_") for p in p_list]
    d_list=sort(unique(df[!,"d"]))
    for d in d_list
        y_err=[df[j,"sem_mi"] for p in p_list for j in 1:length(df[!,"sem_mi"]) if df[j,"d"] == d && df[j,"p"]==p]
        plot!(p_list,[df[j,"avg_mi"] for p in p_list for j in 1:length(df[!,"avg_mi"]) if df[j,"d"] == d && df[j,"p"]==p],yerr=y_err,label = "d=$d,q=$q",xlabel = L"p_{in}",ylabel=L"\langle I(A:B) \rangle / n",marker=(:circle,4), grid=false, legend=:bottomleft,  xtickfontsize=15,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,legendfontsize=14, linewidth=2,xticks=([0.01,0.02,0.03],[0.01,0.02,0.03]))
    end
    savefig("plots/mi_L_$(l)_r_$(r).pdf")
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--size", "-l"
            arg_type = Int64
            required = true
        "--rate", "-r"
            arg_type = Int64
            required = true
    end
    return parse_args(s)
end
if isinteractive() == false
    plot_mi_w_errbar(parse_commandline())
end
