using ArgParse, CSV, DataFrames, JuMP, Ipopt, Plots, LaTeXStrings

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

#given list of layers and starting state produce code
function gen_stabilizers(m,layer_list)
    for layer in layer_list
        for (gate,loc) in layer
            apply!(m,gate,loc)
        end
    end
    return m
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
    #return MixedDestabilizer(stabilizerview(m)[:,1:length(locs)])
end
function erasure_errors(m,p,L,indices)
    if length(stabilizerview(m)) == 0
        return one(MixedDestabilizer,0,L)
    end
    if stabilizerview(m)[1].nqubits!=L
        throw(DomainError(x, "inconsistent L and m"))
    end
    swap = C"_X X_ _Z Z_"
    #n = stabilizerview(m)[1].nqubits
    n = length(indices)
    d = Bernoulli(p)
    locs = rand(d, n)
    j=L
    for (i, s) in enumerate(locs)
        if s == 1
            e = one(MixedDestabilizer,0,1)
            m = m ⊗ e
            j+=1
            apply!(m,swap,[indices[i],j], phases=false)
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
    #return MixedDestabilizer(stabilizerview(m)[:,1:length(locs)])
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
    #println("after bit flip checks: $(stabilizerview(m))")
    #check phase flip errors

    #println("after phase flip checks: $(stabilizerview(m))")
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

function gen_css_code(L, d)
    gateset = gen_css_gates(2)
    m = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    stabilizer_ind = [i for i in 1:L if i%3 != 0]
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        while i < L
            apply!(m, sample(gateset), [i, i+1], phases=false)
            i += 2
        end
        if i == L
            apply!(m, sample(gateset), [L, 1], phases=false)
        end
    end
    return m, stabilizer_ind
end
function distill_0_alt_v7(layer_list,q,p,L,order,input)
    # base case
    if q == 0
        return copy(input)
    else
        if (q-1)==0
            m1 = distill_0_alt_v7(layer_list,q-1,p,L,order,input)
            m2 = distill_0_alt_v7(layer_list,q-1,p,L,order,input)
            m1 = mid_circuit_errors(m1,layer_list,p,L)
            m2 = mid_circuit_errors(m2,layer_list,p,L)
        else
            m1 = distill_0_alt_v7(layer_list,q-1,p,L,order,input)
            m2 = distill_0_alt_v7(layer_list,q-1,p,L,order,input)
        end
    end
    m=steane_verification_single_alt_noisy(m1,m2,p,L,q%2==order)
    return m
end
function gen_bell_pair3(layer_list,L,m0_plus,m0_zero)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    input_plus = gen_stabilizers(copy(m0_plus),layer_list)
    # println("input_plus")
    # println(stabilizerview(input_plus))
    #zero state
    input_zero = gen_stabilizers(copy(m0_zero),layer_list)
    # println("input_zero")
    # println(stabilizerview(input_zero))
    input = input_plus ⊗ input_zero

    for i in 1:L
        apply!(input, CNOT, [i,i+L])
    end

    return input
end
function ec_gadget(m,p,L)
    CNOT = C"XX _X Z_ ZZ"
    # check bit flip errors
    for i in (L+1):(2*L)
        apply!(m,CNOT,[i,i+L])
    end
    m=erasure_errors(m,p,4*L,[i for i in ((2*L)+1):(3*L)])
    for i in ((2*L)+1):(3*L)
        project!(m,single_z(4*L,i),phases=false)
    end
    # check phase flip errors
    for i in (L+1):(2*L)
        apply!(m,CNOT,[(2*L)+i,i])
    end
    m=erasure_errors(m,p,4*L,[i for i in ((3*L)+1):(4*L)])
    for i in ((3*L)+1):(4*L)
        project!(m,single_x(4*L,i),phases=false)
    end

    traceout!(m,[i for i in ((2*L)+1):(4*L)])
    #println("after bit flip checks: $(stabilizerview(m))")
    #check phase flip errors

    #println("after phase flip checks: $(stabilizerview(m))")
    if length(stabilizerview(m)) == 0
        return one(MixedDestabilizer,0,2*L)
    end
    s = stabilizerview(m)[:,1:(2*L)]

    canonicalize!(s)
    stabilizers = [i for i in s if (i != zero(PauliOperator, 2*L)) && (i != -1*zero(PauliOperator, 2*L))]
    if length(stabilizers) == 0
        return one(MixedDestabilizer,0,2*L)
    else
        return MixedDestabilizer(Stabilizer(stabilizers))
    end
end
function measure_stabilizers_bell(bell_pair,code,L)
    for s in stabilizerview(code)
        project!(bell_pair,zero(PauliOperator,L) ⊗ s,phases=false)
        #project!(bell_pair,s ⊗ zero(PauliOperator,L) ,phases=false)
    end
    return bell_pair
end
function calc_mi(bell_pair,code,L)
    bell_pair = measure_stabilizers_bell(bell_pair,code,L)
    total_entropy = 2*L - rank(bell_pair)
    alice = copy(bell_pair)
    traceout!(alice,[i for i in (L+1):(2*L)])
    alice_entropy = L-rank(alice)
    bob = copy(bell_pair)
    traceout!(bob,[i for i in 1:L])
    bob_entropy = L-rank(bob)
    mi = (alice_entropy+bob_entropy-total_entropy)/L
    return mi
end
function gen_input(L,rate)
    x_ones_plus = []
    x_ones_zero = []
    j=1
    x=true
    while j <= L
        if j % rate == 0
            push!(x_ones_plus,true)
            push!(x_ones_zero,false)
        else
            push!(x_ones_plus,x)
            push!(x_ones_zero,x)
            x = !x
        end
        j+=1
    end
    m0_plus= MixedDestabilizer(Stabilizer([ifelse(x,single_x(L,i),single_z(L,i)) for (i,x) in enumerate(x_ones_plus)]))
    m0_zero= MixedDestabilizer(Stabilizer([ifelse(x,single_x(L,i),single_z(L,i)) for (i,x) in enumerate(x_ones_zero)]))
    m0_stabilizers = MixedDestabilizer(Stabilizer([ifelse(x,single_x(L,i),single_z(L,i)) for (i,x) in enumerate(x_ones_plus) if i%rate!=0]))
    return m0_plus,m0_zero,m0_stabilizers
end

function main_mi_test5(args)
    L=args["size"]
    rate=args["rate"]
    p = args["error_rate"]
    # filename = args["filename"]
    println("L: $L")
    println("rate: $rate")
    println("p: $p")
    d = args["depth"]
    q = args["levels"]
    println("d: $d, q:$q")
    rounds = 10
    ns = args["num_samples"]
    #order=1 is bit flip first, order=0 is phase flip first
    plus_order = 0
    zero_order = 1
    m0_plus,m0_zero,m0_stabilizers = gen_input(L,rate)
    open("data/mi_L_$(L)_d_$(d)_r_$(rate)_p_$(p)_q_$(q).csv", "a") do file
        if args["header"]
            write(file, "L,d,r,p,q,samples,$(join(["mi_$(i)" for i in 0:rounds],","))\n")
        end
        mi_tot_list = zeros(rounds+1)
        for k in 1:ns
            println("k: $k")
            layer_list = gen_css_code_layers(L, d)
            code = gen_stabilizers(copy(m0_stabilizers),layer_list)
            bell_pair = gen_bell_pair3(layer_list,L,copy(m0_plus),copy(m0_zero))

            mi_list = zeros(rounds+1)
            mi_list[1] = calc_mi(copy(bell_pair),code,L)
            for i in 1:rounds
                m1 = distill_0_alt_v7(layer_list,q,p,L,plus_order,copy(m0_plus))
                m2 = distill_0_alt_v7(layer_list,q,p,L,zero_order,copy(m0_zero))
                bell_pair = bell_pair ⊗ m1 ⊗ m2
                bell_pair = ec_gadget(bell_pair,p,L)

                mi = calc_mi(copy(bell_pair),code,L)
                mi_list[i+1] = mi
            end
            write(file, "$L,$d,$rate,$p,$q,$k,$(join(mi_list,","))\n")
            mi_tot_list = mi_tot_list .+ mi_list
            if k % 100 == 0
                println("avg_mi: $(join(mi_tot_list ./ k,","))")
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
        "--error_rate", "-p"
            arg_type = Float64
            required = true
        "--rate", "-r"
            arg_type = Int64
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
    main_mi_test5(parse_commandline())
end
