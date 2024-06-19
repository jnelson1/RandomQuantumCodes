using QuantumClifford,SimpleGF2,LinearAlgebra,Random,Distributions,Plots,CSV,DataFrames,ArgParse

function all_functions(finite_set, nx, ny)
  I = Iterators.product(fill(finite_set, nx*ny)...)
  return (reshape(collect(i), (nx,ny)) for i in I)
end
function gen_css_gates(dim)
    fs = [0, 1]
    gates = []
    for Uxx in all_functions(fs, dim, dim)
        Uxxgf2 = map(GF2, Uxx)
        if (det(Uxxgf2)).val == 1
            Uxxgf2inv = inv(Uxxgf2)
            Uzzgf2 = transpose(Uxxgf2inv)
            push!(gates, convert(Uxxgf2, Uzzgf2, dim))
        end
    end
    return gates
end
function convert(Uxx, Uzz, n)
    # xcolumns = [PauliOperator(0x0, [Uxx[i, j].val==1 for i in 1:n], [false for i in 1:n]) for j in 1:n]
    # zcolumns = [PauliOperator(0x0, [false for i in 1:n], [Uzz[i, j].val==1 for i in 1:n]) for j in 1:n]
    # return CliffordColumnForm(cat(xcolumns, zcolumns, dims=1))
    Cstring_xlist = [join([ifelse(Uxx[i, j].val==1,"X","_") for i in 1:n]) for j in 1:n]
    Cstring_zlist = [join([ifelse(Uzz[i, j].val==1,"Z","_") for i in 1:n]) for j in 1:n]
    Cstring_list = vcat(Cstring_xlist,Cstring_zlist)
    Cstring = join(Cstring_list," ")
    exp = Meta.parse("gate = C\"$(Cstring)\"")
    eval(exp)
    return gate

end
function gen_css_code_layers(L,d)
    gateset = gen_css_gates(2)
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            push!(layer,(sample(gateset),[i,i+1]))
            i += 2
        end
        if i == L
            push!(layer,(sample(gateset), [L, 1]))
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function erasure_errors(m,p,L)
    if length(stabilizerview(m)) == 0
        return one(MixedDestabilizer,0,L)
    end
    if stabilizerview(m)[1].nqubits!=L
        throw(DomainError(x, "inconsistent L and m"))
    end
    swap = C"_X X_ _Z Z_"
    d = Bernoulli(p)
    locs = rand(d, L)
    j=L
    for (i, s) in enumerate(locs)
        if s == 1
            e = one(MixedDestabilizer,0,1)
            m = m ⊗ e
            j+=1
            apply!(m,swap,[i,j], phases=false)
            traceout!(m,[j])
        end
    end
    if length(stabilizerview(m)) == 0
        return one(MixedDestabilizer,0,L)
    end
    s = stabilizerview(m)[:,1:L]
    canonicalize!(s)
    stabilizers = [i for i in s if (i != zero(PauliOperator, L)) && (i != -1*zero(PauliOperator, L))]
    if length(stabilizers) == 0
        return one(MixedDestabilizer,0,L)
    else
        return MixedDestabilizer(Stabilizer(stabilizers))
    end
end
function mid_circuit_errors(m,layer_list,p,L)
    for layer in layer_list
        for (gate,loc) in layer
            apply!(m,gate,loc)
        end
        m = erasure_errors(m,p,L)
    end
    return m
end
function steane_verification_single_alt_noisy(m1,m2,p,L,alt)
    CNOT = C"XX _X Z_ ZZ"
    m = m1 ⊗ m2
    if alt == true
        # check bit flip errors
        for i in 1:L
            apply!(m,CNOT,[i,i+L])
        end
        m=erasure_errors(m,p,2*L)
        for i in 1:L
            project!(m,single_z(2*L,L+i),phases=false)
        end
    else
    # check phase flip errors
        for i in 1:L
            apply!(m,CNOT,[L+i,i])
        end
        m=erasure_errors(m,p,2*L)
        for i in 1:L
            project!(m,single_x(2*L,L+i),phases=false)
        end
    end
    traceout!(m,[i for i in (L+1):(2*L)])
    #check phase flip errors

    if length(stabilizerview(m)) == 0
        return one(MixedDestabilizer,0,L)
    end
    s = stabilizerview(m)[:,1:L]

    canonicalize!(s)
    stabilizers = [i for i in s if (i != zero(PauliOperator, L)) && (i != -1*zero(PauliOperator, L))]
    if length(stabilizers) == 0
        return one(MixedDestabilizer,0,L)
    else
        return MixedDestabilizer(Stabilizer(stabilizers))
    end
end
#only one bit flip or phase flip is checked at a time
#noisy encoding and noisy CNOTs interleaved
#q is level
function distill_0_alt_v5(layer_list,q,p,L)
    # base case
    entropy_list = []
    if q == 0
        return MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L])),entropy_list
    else
        if (q-1)==0
            m1,entropy_list1 = distill_0_alt_v5(layer_list,q-1,p,L)
            m2,entropy_list2 = distill_0_alt_v5(layer_list,q-1,p,L)
            entropy_list = vcat(entropy_list1,entropy_list2)
            m1 = mid_circuit_errors(m1,layer_list,p,L)
            m2 = mid_circuit_errors(m2,layer_list,p,L)
        else
            m1,entropy_list1 = distill_0_alt_v5(layer_list,q-1,p,L)
            m2,entropy_list2 = distill_0_alt_v5(layer_list,q-1,p,L)
            entropy_list = vcat(entropy_list1,entropy_list2)
        end
    end
    m=steane_verification_single_alt_noisy(m1,m2,p,L,q%2==0)
    entropy = (L - rank(m))/L
    push!(entropy_list, (q,entropy))
    return m,entropy_list
end
# noisy encoding interweaved with noisy CNOTS
function main_0_alt_single(args)
    L = args["size"]
    d=args["depth"]
    p = args["error_rate"]
    num_samples = args["num_samples"]
    levels = args["levels"]
    open("data/$(args["filename"]).csv", "a") do file
        if args["header"] == true
            write(file, "L,d,p,levels,num_samples,$(join(["avg_entropy_$(i)" for i in 1:levels],",")),$(join(["sem_entropy_$(i)" for i in 1:levels],","))\n")
        end
        entropy_dict = Dict([(q,[]) for q in 0:levels])
        for j in 1:num_samples
            println("j: $(j)")
            layer_list = gen_css_code_layers(L, d)
            m,entropy_list = distill_0_alt_v5(layer_list,levels,p,L)
            for (q,e) in entropy_list
                push!(entropy_dict[q],e)
            end
        end
        entropy_avg = [sum(entropy_dict[q])/length(entropy_dict[q]) for q in 1:levels]
        entropy_std = [sqrt(sum((entropy_dict[q] .- entropy_avg[q]).^2)/(length(entropy_dict[q]) - 1)) for q in 1:levels]
        entropy_sem = [entropy_std[q] / sqrt(length(entropy_dict[q])) for q in 1:levels]
        write(file, "$L,$d,$p,$levels,$num_samples,$(join(entropy_avg,",")),$(join(entropy_sem,","))\n")
    end
end
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--filename", "-f"
            required = true
        "--size", "-l"
            arg_type = Int64
            required = true
        "--error_rate", "-p"
            arg_type = Float64
            required = true
        "--num_samples", "-n"
            arg_type = Int64
            required = true
        "--depth", "-d"
            arg_type = Int64
            required = true
        "--levels", "-q"
            arg_type = Int64
            required = true
        "--header", "-a"
            action = :store_true
    end
    return parse_args(s)
end
if isinteractive() == false
    main_0_alt_single(parse_commandline())
end
