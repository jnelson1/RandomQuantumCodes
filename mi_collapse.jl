using ArgParse, CSV, DataFrames, JuMP, Ipopt, Plots, LaTeXStrings

function main_mi(args)
    L=args["size"]
    r=args["rate"]

    df = DataFrame(CSV.File("data/mi_L_$(L)_r_$(r).csv"))
    m = size(df)[1]
    d_list = []
    pin_list = []
    pout_list = []
    sem_list = []
    for i in 1:m
        if df[i,"p"] >=0.015 && df[i,"p"] <=0.025
            push!(d_list,df[i,"d"])
            push!(pin_list,df[i,"p"])
            push!(pout_list,df[i,"avg_mi"])
            push!(sem_list,df[i,"sem_mi"])
        end
    end
    data = hcat(d_list,pin_list,pout_list,sem_list)
    n = size(data)[1]
    p_c_list = zeros(Float64, n)
    λ_list = zeros(Float64, n)
    A_list = zeros(Float64, n)
    B_list = zeros(Float64, n)
    C_list = zeros(Float64, n)
    for i in 1:n
        new_data = data[1:n .!= i,:]
        p_c, λ, A, B, C = fit(new_data, 0.0)
        p_c_list[i] = p_c
        λ_list[i] = λ
        A_list[i] = A
        B_list[i] = B
        C_list[i] = C
    end
    p_c_avg = sum(p_c_list)/n
    λ_avg = sum(λ_list)/n
    A_avg = sum(A_list)/n
    B_avg = sum(B_list)/n
    C_avg = sum(C_list)/n
    variance = ((n-1)/n)*sum((p_c_list .- p_c_avg).^2)
    std = sqrt(variance)
    #p_c, λ, A, B = fit(data, pc_guess)
    println("p_c_avg: $(p_c_avg)")
    println("λ_avg: $(λ_avg)")
    println("A_avg: $(A_avg)")
    println("B_avg: $(B_avg)")
    println("C_avg: $(C_avg)")
    println("std: $(std)")
    #my_plot(p_c_avg, λ_avg, data, args["fileid"])
    my_plot_collapse_mi(p_c_avg, λ_avg, data, "mi_L_$(L)_r_$(r)")
    #my_plot_cross(data, args["fileid"])
end
function fit(data, pc_guess)
    n = size(data)[1]
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    @variable(model, p_c, start=pc_guess)
    @variable(model, λ)
    @variable(model, A)
    @variable(model, B)
    @variable(model, C)
    @NLexpression(model, mse, sum((A+(B*(data[i,2]-p_c)*(data[i,1]^λ))+(C*((data[i,2]-p_c)*(data[i,1]^λ))^2) - data[i,3])^2 for i in 1:n))
    #@NLexpression(model, mse, sum((A+(B*(data[i,2]-p_c)*(data[i,1]^λ)) - data[i,3])^2 for i in 1:n))
    #@NLexpression(model, mse, sum((A+(B*(data[i,2]-p_c)*(exp(data[i,1]*λ)))+(C*((data[i,2]-p_c)*(exp(data[i,1]*λ)))^2) - data[i,3])^2 for i in 1:n))
    #@NLexpression(model, mse, sum((A+(B*(data[i,2]-p_c)*(data[i,1]^λ)) - data[i,3])^2 for i in 1:n))
    @NLobjective(model, Min, mse)
    optimize!(model)
    return value(p_c), value(λ), value(A), value(B), value(C)
    #return value(p_c), value(λ), value(A), value(B)
end
function my_plot_collapse_mi(p_c, λ, data, fileid)
    plot()
    n = size(data)[1]
    dlist = sort(unique([data[i,1] for i in 1:n]))
    # r=4
    # for d in [2,4,6]
    #     p_in = [data[i,2] for i in 1:n if data[i,1]==d]
    #     p_out = [data[i,3] for i in 1:n if data[i,1]==d]
    #     sem=[data[i,4] for i in 1:n if data[i,1]==d]
    #     x = (p_in .- p_c) .* d^λ
    #     scatter!(x, p_out, yerr = 2*sem, label="d=$(d),q=$(d)",marker=(:circle,4), grid=false, legend=:bottomleft, xlabel = L"(p_{in}-p_c)d^\lambda", ylabel =  L"I(A:B)", xtickfontsize=15,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,legendfontsize=14, linewidth=2,margin=5Plots.mm,yticks=([0.20,0.25,0.3,0.35],[0.20,0.25,0.3,0.35]))
    #     #xticks=([0.0,0.02],[0.0,0.02])
    #     #yticks=([0.35,0.4],[0.35,0.4])
    # end
    # r=5
    # for d in [2,4,6]
    #     p_in = [data[i,2] for i in 1:n if data[i,1]==d]
    #     p_out = [data[i,3] for i in 1:n if data[i,1]==d]
    #     sem=[data[i,4] for i in 1:n if data[i,1]==d]
    #     x = (p_in .- p_c) .* d^λ
    #     scatter!(x, p_out, yerr = 2*sem, label="d=$(d),q=$(d)",marker=(:circle,4), grid=false, legend=:bottomleft, xlabel = L"(p_{in}-p_c)d^\lambda", ylabel =  L"I(A:B)", xtickfontsize=15,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,legendfontsize=14, linewidth=2,margin=5Plots.mm,yticks=([0.2,0.25],[0.2,0.25]))
    #     #ylims=(0.275,0.35),yticks=([0.3,0.35],[0.3,0.35])
    # end
    # r=10
    # for d in [2,4,6]
    #     p_in = [data[i,2] for i in 1:n if data[i,1]==d]
    #     p_out = [data[i,3] for i in 1:n if data[i,1]==d]
    #     sem=[data[i,4] for i in 1:n if data[i,1]==d]
    #     x = (p_in .- p_c) .* d^λ
    #     scatter!(x, p_out, yerr = 2*sem, label="d=$(d),q=$(d)",marker=(:circle,4), grid=false, legend=:bottomleft, xlabel = L"(p_{in}-p_c)d^\lambda", ylabel =  L"I(A:B)", xtickfontsize=15,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,legendfontsize=14, linewidth=2,margin=5Plots.mm,ylims=(0.05,0.15),yticks=([0.05,0.1,0.15],[0.05,0.1,0.15]))
    # end
    # r=3
    for d in dlist
        p_in = [data[i,2] for i in 1:n if data[i,1]==d]
        p_out = [data[i,3] for i in 1:n if data[i,1]==d]
        sem=[data[i,4] for i in 1:n if data[i,1]==d]
        x = (p_in .- p_c) .* d^λ
        scatter!(x, p_out, yerr = 2*sem, label="d=$(d),q=$(d)",marker=(:circle,4), grid=false, legend=false, xlabel = L"(p_{in}-p_c)d^\lambda", xtickfontsize=22,ytickfontsize=22,xguidefontsize=24, linewidth=2,margin=5Plots.mm,xticks=([-0.01,0.0,0.01,0.02],[-0.01,0.0,0.01,0.02]),ylims=(0.29,0.5),yticks=([0.3,0.4,0.5],[0.3,0.4,0.5]))
    end
    # for inset
    # for d in [2,4,6]
    #     p_in = [data[i,2] for i in 1:n if data[i,1]==d]
    #     p_out = [data[i,3] for i in 1:n if data[i,1]==d]
    #     sem=[data[i,4] for i in 1:n if data[i,1]==d]
    #     x = (p_in .- p_c) .* d^λ
    #     scatter!(x, p_out, yerr = 2*sem, label="d=$(d),q=$(d)",marker=(:circle,4), grid=false, legend=false, xlabel = L"(p_{in}-p_c)d^\lambda", ylabel =  L"I(A:B)", xtickfontsize=22,ytickfontsize=22,xguidefontsize=24,yguidefontsize=24, linewidth=2,margin=5Plots.mm,xticks=([-0.01,0.0,0.01,0.02],[-0.01,0.0,0.01,0.02]),ylims=(0.29,0.5),yticks=([0.3,0.4,0.5],[0.3,0.4,0.5]))
    #     #yticks=([0.5,0.55],[0.5,0.55])
    # end

    # r=50
    # for d in [2,4,6]
    #     p_in = [data[i,2] for i in 1:n if data[i,1]==d]
    #     p_out = [data[i,3] for i in 1:n if data[i,1]==d]
    #     sem=[data[i,4] for i in 1:n if data[i,1]==d]
    #     x = (p_in .- p_c) .* d^λ
    #     scatter!(x, p_out, yerr = 2*sem, label="d=$(d),q=$(d)",marker=(:circle,4), grid=false, legend=:bottomleft, xlabel = L"(p_{in}-p_c)d^\lambda", ylabel =  L"I(A:B)", xtickfontsize=15,ytickfontsize=16,xguidefontsize=18,yguidefontsize=18,legendfontsize=14, linewidth=2,margin=5Plots.mm,xticks=([-0.01,0.0,0.01],[-0.01,0.0,0.01]),yticks=([0.02,0.025],[0.02,0.025]))
    #     #ylims=(0.0275,0.035),xticks=([-0.01,0.0,0.01,0.02],[-0.01,0.0,0.01]),yticks=([0.03,0.035],[0.03,0.035])
    # end
    savefig("plots/thresh_est_$(fileid)_v2.pdf")

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
    main_mi(parse_commandline())
end
