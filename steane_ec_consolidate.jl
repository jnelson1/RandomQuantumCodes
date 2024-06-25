using CSV,DataFrames,ArgParse,Statistics

function consolidate(args)
    L=args["size"]
    r=args["rate"]
    dlist = args["dlist"]
    qlist = args["qlist"]
    plist = args["plist"]

    open("data/mi_L_$(L)_r_$(r).csv", "a") do file
        write(file, "L,d,r,p,num_samples,avg_mi,sem_mi\n")
        for (d,q) in zip(dlist,qlist)
            for p in plist
                filename = "data/mi_L_$(L)_d_$(d)_r_$(r)_p_$(p)_q_$(q).csv"
                df = CSV.read(filename, DataFrame)
                ns=length(df[!,"samples"])
                avg_mi = mean(df[!,"mi_10"])
                sem_mi = (std(df[!,"mi_10"]))/sqrt(ns)
                write(file,"$L,$d,$r,$p,$ns,$avg_mi,$sem_mi\n")
            end
        end
    end
end
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--size", "-l"
            arg_type = Int64
            required = true
        "--plist", "-p"
            arg_type = Float64
            required = true
            nargs='*'
        "--rate", "-r"
            arg_type = Int64
            required = true
        "--dlist", "-d"
            arg_type = Int64
            required = true
            nargs='*'
        "--qlist", "-q"
            arg_type = Int64
            required = true
            nargs='*'
end
    return parse_args(s)
end
if isinteractive() == false
    consolidate(parse_commandline())
end
