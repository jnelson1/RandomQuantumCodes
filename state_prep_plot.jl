using Plots,CSV,DataFrames,ArgParse,LaTeXStrings

function plot_entropy6(args)
    filename = "data/$(args["filename"]).csv"
    println(filename)
    df = CSV.read(filename, DataFrame)
    p_list = sort(unique(df[!,"p"]))
    rounds = args["rounds"]
    plot()
    for i in rounds
        y_err=[df[j,"sem_entropy_$i"] for p in p_list for j in 1:length(df[!,"sem_entropy_$i"]) if df[j,"p"]==p]
        plot!(p_list,[df[j,"avg_entropy_$i"] for p in p_list for j in 1:length(df[!,"avg_entropy_$i"]) if df[j,"p"]==p],yerr=y_err,label = "q=$i",xlabel = L"p_{in}",ylabel=L"\langle S(\rho) \rangle /n",marker=(:circle,4), grid=false, legend=:topleft,  xtickfontsize=15,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,legendfontsize=14, linewidth=2,xticks=([0.01,0.02,0.03],[0.01,0.02,0.03]))
        #xticks=([0.01,0.02,0.03],[0.01,0.02,0.03]),yticks=([0.025,0.03,0.035],[0.025,0.03,0.035])
    end
    savefig("plots/$(args["filename"]).pdf")
end
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--filename", "-f"
            required = true
        "--rounds", "-r"
            arg_type = Int64
            nargs='*'
            required=true
    end
    return parse_args(s)
end

if isinteractive() == false
    plot_entropy6(parse_commandline())
end
