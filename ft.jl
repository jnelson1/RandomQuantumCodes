using QuantumClifford,SimpleGF2,LinearAlgebra,Random,Distributions,Plots,CSV,DataFrames

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
function gen_noncss_code_v2_layers(L,d)
    CNOT = C"XX _X Z_ ZZ"
    cnot_gateset = [CNOT, permute(CNOT, [2,1])]
    Hadamard = C"Z X"
    layer_list = []
    for l in 0:(d - 1)
        layer = []
        i = 1 + (l % 2)
        while i < L
            push!(layer,(sample(cnot_gateset), [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer,(sample(cnot_gateset), [L, 1]))
        end
        for j in 1:L
            if sample([0,1]) == 1
                push!(layer,(Hadamard, [j]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_css_code_v3_layers(L,d)
    CNOT = C"XX _X Z_ ZZ"
    cnot_gateset = [CNOT, permute(CNOT, [2,1])]
    layer_list = []
    for l in 0:(d - 1)
        layer = []
        i = 1 + (l % 2)
        while i < L
            push!(layer,(sample(cnot_gateset), [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer,(sample(cnot_gateset), [L, 1]))
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_noncss_code_v4_layers(L, d)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            push!(layer, (CNOT, [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer, (CNOT, [L, 1]))
        end
        for j in 1:L
            if sample([0,1]) == 1
                push!(layer,(Hadamard, [j]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_noncss_code_v5_layers(L, d, p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    coin = Bernoulli(p)
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            push!(layer, (CNOT, [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer, (CNOT, [L, 1]))
        end
        for j in 1:L
            if rand(coin) == 1
                push!(layer,(Hadamard, [j]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_noncss_code_v5(L, d, p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    coin = Bernoulli(p)
    m = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        while i < L
            apply!(m, CNOT, [i, i+1])
            i += 2
        end
        if i == L
            apply!(m, CNOT, [L, 1])
        end
        for j in 1:L
            if rand(coin) == 1
                apply!(m, Hadamard, [j])
            end
        end
    end
    return m
end
function gen_css_code_v5_layers(L, d, p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    coin = Bernoulli(p)
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            push!(layer, (CNOT, [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer, (CNOT, [L, 1]))
        end
        if rand(coin) == 1
            for j in 1:L
                push!(layer,(Hadamard, [j]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_css_code_v6_layers(L, d)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            push!(layer, (CNOT, [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer, (CNOT, [L, 1]))
        end
        if l == 0
            for j in 1:L
                push!(layer,(Hadamard, [j]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_css_code_v5(L, d, p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    coin = Bernoulli(p)
    m = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        while i < L
            apply!(m, CNOT, [i, i+1])
            i += 2
        end
        if i == L
            apply!(m, CNOT, [L, 1])
            #push!(layer, (CNOT, [L, 1]))
        end
        if rand(coin) == 1
            for j in 1:L
                apply!(m, Hadamard, [j])
            end
        end
    end
    return m
end
function gen_css_code_v7_layers(L, d, p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    gateset = gen_css_gates(2)
    coin = Bernoulli(p)
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            if rand(coin) == 1
                push!(layer,(sample(gateset), [i, i+1]))
            else
                push!(layer, (CNOT, [i, i+1]))
            end
            i += 2
        end
        if i == L
            if rand(coin)== 1
                push!(layer,(sample(gateset), [L, 1]))
            else
                push!(layer, (CNOT, [L, 1]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_css_code_v8_layers(L, d)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    gateset = gen_css_gates(2)
    #coin = Bernoulli(p)
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            if l == 0
                push!(layer,(sample(gateset), [i, i+1]))
            else
                push!(layer, (CNOT, [i, i+1]))
            end
            i += 2
        end
        if i == L
            if l== 0
                push!(layer,(sample(gateset), [L, 1]))
            else
                push!(layer, (CNOT, [L, 1]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_css_code_v9_layers(L, d, p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    gateset = gen_css_gates(2)
    coin = Bernoulli(p)
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            if l == 0 && rand(coin) == 1
                push!(layer,(sample(gateset), [i, i+1]))
            else
                push!(layer, (CNOT, [i, i+1]))
            end
            i += 2
        end
        if i == L
            if l== 0 && rand(coin) == 1
                push!(layer,(sample(gateset), [L, 1]))
            else
                push!(layer, (CNOT, [L, 1]))
            end
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_css_code_v10_layers(L, d, p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    Swap = C"_X X_ _Z Z_"
    gateset = gen_css_gates(2)
    coin = Bernoulli(p)
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            if l == 0 && rand(coin) == 1
                push!(layer,(Swap, [i, i+1]))
            end
            push!(layer, (CNOT, [i, i+1]))
            i += 2
        end
        if i == L
            if l== 0 && rand(coin) == 1
                push!(layer,(Swap, [L, 1]))
            end
            push!(layer, (CNOT, [L, 1]))
        end
        push!(layer_list,layer)
    end
    return layer_list
end
# function gen_css_code_v2(L, d)
#     CNOT = C"XX _X Z_ ZZ"
#     cnot_gateset = [CNOT, permute(CNOT, [2,1])]
#     Hadamard = C"Z X"
#     m = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
#     stabilizer_ind = [i for i in 1:L if i%3 != 0]
#     for l in 0:(d - 1)
#         i = 1 + (l % 2)
#         while i < L
#             apply!(m, sample(cnot_gateset), [i, i+1], phases=false)
#             i += 2
#         end
#         if i == L
#             apply!(m, sample(cnot_gateset), [L, 1], phases=false)
#         end
#         for j in 1:L
#             if sample([0,1]) == 1
#                 apply!(m, Hadamard, [j])
#             end
#         end
#     end
#     return m, stabilizer_ind
# end
function gen_noncss_code_layers(L, d)
    layer_list = []
    for l in 0:(d - 1)
        layer = []
        i = 1 + (l % 2)
        while i < L
            push!(layer,(random_clifford(2), [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer,(random_clifford(2), [L, 1]))
        end
        push!(layer_list,layer)
    end
    return layer_list
end
# function gen_noncss_code(L, d)
#     m = MixedDestabilizer(Stabilizer([single_z(L, i) for i in 1:L]))
#     stabilizer_ind = [i for i in 1:L if i%3 != 0]
#     for l in 0:(d - 1)
#         i = 1 + (l % 2)
#         while i < L
#             apply!(m, random_clifford(2), [i, i+1])
#             i += 2
#         end
#         if i == L
#             apply!(m, random_clifford(2), [L, 1])
#         end
#     end
#     return m, stabilizer_ind
# end
function gen_cliff_east_code_layers(L, d)
    CNOT = C"XX _X Z_ ZZ"
    layer_list = []
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        layer = []
        while i < L
            push!(layer, (CNOT, [i, i+1]))
            i += 2
        end
        if i == L
            push!(layer, (CNOT, [L, 1]))
        end
        push!(layer_list,layer)
    end
    return layer_list
end
function gen_cliff_east_code(L, d)
    CNOT = C"XX _X Z_ ZZ"
    m = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    stabilizer_ind = [i for i in 1:L if i%3 != 0]
    for l in 0:(d - 1)
        i = 1 + (l % 2)
        while i < L
            apply!(m, CNOT, [i, i+1])
            i += 2
        end
        if i == L
            apply!(m, CNOT, [L, 1])
        end
    end
    return m, stabilizer_ind
end
#assume block_size is even
#assume L is 2^k * block_size
#assume input qubits are repeating pattern of stabilizer,logical,stabilizer
function gen_cliff_east_code_periodic_block_layers(L,t,block_size)
    CNOT = C"XX _X Z_ ZZ"
    layer_list = []
    while block_size <= L
        for l in 0:(t-1)
            i = 1 + (l % 2)
            layer = []
            while i <= L
                if i % block_size == 0
                    push!(layer, (CNOT, [i, i-block_size+1]))
                    i+=2
                else
                    push!(layer, (CNOT, [i, i+1]))
                    i += 2
                end
            end
            push!(layer_list,layer)
        end
        block_size = 2*block_size
    end
    return layer_list
end
function gen_cliff_east_code_periodic_block_code(L,t,block_size)
    CNOT = C"XX _X Z_ ZZ"
    m = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    stabilizer_ind = [i for i in 1:L if (i+1)%3 != 0]
    while block_size <= L
        for l in 0:(t-1)
            i = 1 + (l % 2)
            while i <= L
                if i % block_size == 0
                    apply!(m,CNOT,[i, i-block_size+1])
                    i+=2
                else
                    apply!(m,CNOT,[i, i+1])
                    i+=2
                end
            end
        end
        block_size = 2*block_size
    end
    return m,stabilizer_ind
end
function gen_cliff_east_code_open_block_layers(num_blocks,t,block_size)
    CNOT = C"XX _X Z_ ZZ"
    L=num_blocks*block_size
    layer_list = []
    while block_size <= L
        for l in 0:(t-1)
            i = 1 + (l % 2)
            layer = []
            block_end_ind = block_size
            while block_end_ind <= L
                while i < block_end_ind
                    push!(layer, (CNOT, [i, i+1]))
                    i += 2
                end
                i = 1+(l % 2)+block_end_ind
                block_end_ind += block_size
            end
            push!(layer_list,layer)
        end
        block_size = 2*block_size
    end
    return layer_list
end
function gen_cliff_east_code_open_block_code(num_blocks,t,block_size)
    CNOT = C"XX _X Z_ ZZ"
    L=num_blocks*block_size
    stabilizers = []
    midpoint = (block_size ÷ 2)+1
    shift = 0
    stabilizer_ind = []
    for _ in 1:num_blocks
        for i in 1:(midpoint-1)
            push!(stabilizer_ind,i+shift)
            if i % 2 == 1
                push!(stabilizers,single_x(L,i+shift))
            else
                push!(stabilizers,single_z(L,i+shift))
            end
        end
        push!(stabilizers, single_z(L,midpoint+shift))
        for i in (midpoint+1):block_size
            push!(stabilizer_ind,i+shift)
            if i % 2 == 1
                push!(stabilizers,single_z(L,i+shift))
            else
                push!(stabilizers,single_x(L,i+shift))
            end
        end
        shift+=block_size
    end
    m = MixedDestabilizer(Stabilizer([s for s in stabilizers]))
    while block_size <= L
        #println("block_size: $block_size")
        for l in 0:(t-1)
            #println("l: $l")
            i = 1 + (l % 2)
            block_end_ind = block_size
            while block_end_ind <= L
                while i < block_end_ind
                    apply!(m,CNOT,[i, i+1])
                    i += 2
                end
                i = 1+(l % 2)+block_end_ind
                block_end_ind += block_size
            end
            #println("m: $(stabilizerview(m))")
        end
        block_size = 2*block_size
    end
    return m,stabilizer_ind
end
function gen_randcss_code_open_block_layers(num_blocks,t,block_size)
    gateset = gen_css_gates(2)
    L=num_blocks*block_size
    layer_list = []
    while block_size <= L
        for l in 0:(t-1)
            i = 1 + (l % 2)
            layer = []
            block_end_ind = block_size
            while block_end_ind <= L
                while i < block_end_ind
                    push!(layer, (sample(gateset), [i, i+1]))
                    i += 2
                end
                i = 1+(l % 2)+block_end_ind
                block_end_ind += block_size
            end
            push!(layer_list,layer)
        end
        block_size = 2*block_size
    end
    return layer_list
end
function gen_randcss_code_open_block_code(num_blocks,t,block_size)
    gateset = gen_css_gates(2)
    L=num_blocks*block_size
    stabilizers = []
    midpoint = (block_size ÷ 2)+1
    shift = 0
    stabilizer_ind = []
    for _ in 1:num_blocks
        for i in 1:(midpoint-1)
            push!(stabilizer_ind,i+shift)
            if i % 2 == 1
                push!(stabilizers,single_x(L,i+shift))
            else
                push!(stabilizers,single_z(L,i+shift))
            end
        end
        push!(stabilizers, single_z(L,midpoint+shift))
        for i in (midpoint+1):block_size
            push!(stabilizer_ind,i+shift)
            if i % 2 == 1
                push!(stabilizers,single_z(L,i+shift))
            else
                push!(stabilizers,single_x(L,i+shift))
            end
        end
        shift+=block_size
    end
    m = MixedDestabilizer(Stabilizer([s for s in stabilizers]))
    while block_size <= L
        #println("block_size: $block_size")
        for l in 0:(t-1)
            #println("l: $l")
            i = 1 + (l % 2)
            block_end_ind = block_size
            while block_end_ind <= L
                while i < block_end_ind
                    apply!(m,sample(gateset),[i, i+1])
                    i += 2
                end
                i = 1+(l % 2)+block_end_ind
                block_end_ind += block_size
            end
            #println("m: $(stabilizerview(m))")
        end
        block_size = 2*block_size
    end
    return m,stabilizer_ind
end
function gen_randcss2_code_open_block_layers(num_blocks,t,block_size,p)
    CNOT = C"XX _X Z_ ZZ"
    Hadamard = C"Z X"
    cnot_gateset = [CNOT, permute(CNOT, [2,1])]
    coin = Bernoulli(p)
    L=num_blocks*block_size
    layer_list = []
    while block_size <= L
        for l in 0:(t-1)
            i = 1 + (l % 2)
            layer = []
            block_end_ind = block_size
            while block_end_ind <= L
                while i < block_end_ind
                    push!(layer, (sample(cnot_gateset), [i, i+1]))
                    i += 2
                end
                i = 1+(l % 2)+block_end_ind
                block_end_ind += block_size
            end
            if rand(coin) == 1
                for j in 1:L
                    push!(layer,(Hadamard, [j]))
                end
            end
            push!(layer_list,layer)
        end
        block_size = 2*block_size
    end
    return layer_list
end
function gen_stabilizers(m,layer_list)
    for layer in layer_list
        for (gate,loc) in layer
            apply!(m,gate,loc)
        end
    end
    return m
end
function mid_circuit_errors(m,layer_list,p)
    for layer in layer_list
        for (gate,loc) in layer
            apply!(m,gate,loc)
        end
        m = erasure_errors(m,p)
    end
    return m
end
function erasure_errors(m,p)
    swap = C"_X X_ _Z Z_"
    n = stabilizerview(m)[1].nqubits
    d = Bernoulli(p)
    locs = rand(d, n)
    for (i, s) in enumerate(locs)
        if s == 1
            e = one(MixedDestabilizer,0,1)
            m = m ⊗ e
            n+=1
            apply!(m,swap,[i,n], phases=false)
            traceout!(m,[n])
        end
    end
    return MixedDestabilizer(stabilizerview(m)[:,1:length(locs)])
end
# function erasure_errors(m,p)
#     swap = C"_X X_ _Z Z_"
#     n = size(m)[2]
#     d = Bernoulli(p)
#     locs = rand(d, n)
#     for (i, s) in enumerate(locs)
#         if s == 1
#             e = zero(Stabilizer, 1)
#             m = m ⊗ e
#             n+=1
#             apply!(m,swap,[i,n], phases=false)
#             traceout!(m,[n])
#         end
#     end
#     m = m[:,1:length(locs)]
#     return m
# end

function measure_stabilizers(code,m,stabilizer_ind)
    # add ancillas to make them the same size
    # n1 = stabilizerview(m)[1].nqubits
    # e = one(MixedDestabilizer, 0, 1)
    # n2 = stabilizerview(code)[1].nqubits
    # while n2 < n1
    #     code = code ⊗ e
    #     n2 += 1
    # end

    for (i,s) in enumerate(stabilizerview(code))
        if i in stabilizer_ind
            project!(m,s,phases=false)
        end
    end
end
function distill(layer_list,L,d,p)
    CNOT = C"XX _X Z_ ZZ"
    # base case
    if d == 0
        return MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    end

    m1 = distill(layer_list[1:(d-1)],L,d-1,p)
    m1 = mid_circuit_errors(m1,[layer_list[d]],L,p)

    m2 = distill(layer_list[1:(d-1)],L,d-1,p)
    m2 = mid_circuit_errors(m2,[layer_list[d]],L,p)

    m3 = distill(layer_list[1:(d-1)],L,d-1,p)
    m3 = mid_circuit_errors(m3,[layer_list[d]],L,p)

    m = m1 ⊗ m2 ⊗ m3
    #println(rank(m))
    # check bit flip errors
    for i in 1:L
        apply!(m,CNOT,[i,i+L])
    end
    for i in 1:L
        project!(m,single_z(3*L,L+i),phases=false)
    end

    #check phase flip errors
    for i in 1:L
        apply!(m,CNOT,[2*L+i,i])
    end
    for i in 1:L
        project!(m,single_x(3*L,2*L+i),phases=false)
    end
    m = MixedDestabilizer(stabilizerview(m)[:,1:L])
    return m
end
function distill2(layer_list,L,q,p)
    CNOT = C"XX _X Z_ ZZ"
    # base case
    if q == 0
        return MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    end

    m1 = distill2(layer_list,L,q-1,p)
    m1 = gen_stabilizers(m1,layer_list[q],L)
    m1 = erasure_errors(m1,p)

    m2 = distill2(layer_list,L,q-1,p)
    m2 = gen_stabilizers(m2,layer_list[q],L)
    m2 = erasure_errors(m2,p)

    m3 = distill2(layer_list,L,q-1,p)
    m3 = gen_stabilizers(m3,layer_list[q],L)
    m3 = erasure_errors(m3,p)

    m = m1 ⊗ m2 ⊗ m3
    #println(rank(m))
    # check bit flip errors
    for i in 1:L
        apply!(m,CNOT,[i,i+L])
    end
    for i in 1:L
        project!(m,single_z(3*L,L+i),phases=false)
    end

    #check phase flip errors
    for i in 1:L
        apply!(m,CNOT,[2*L+i,i])
    end
    for i in 1:L
        project!(m,single_x(3*L,2*L+i),phases=false)
    end
    m = MixedDestabilizer(stabilizerview(m)[:,1:L])
    return m
end

function distill3(layer_list,L,q,p)
    CNOT = C"XX _X Z_ ZZ"
    # base case
    if q == 0
        return MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    end

    m1 = distill3(layer_list,L,q-1,p)
    m1 = mid_circuit_errors(m1,layer_list[q],L,p)

    m2 = distill3(layer_list,L,q-1,p)
    m2 = mid_circuit_errors(m2,layer_list[q],L,p)

    m3 = distill3(layer_list,L,q-1,p)
    m3 = mid_circuit_errors(m3,layer_list[q],L,p)

    m = m1 ⊗ m2 ⊗ m3
    #println(rank(m))
    # check bit flip errors
    for i in 1:L
        apply!(m,CNOT,[i,i+L])
    end
    for i in 1:L
        project!(m,single_z(3*L,L+i),phases=false)
    end

    #check phase flip errors
    for i in 1:L
        apply!(m,CNOT,[2*L+i,i])
    end
    for i in 1:L
        project!(m,single_x(3*L,2*L+i),phases=false)
    end
    m = MixedDestabilizer(stabilizerview(m)[:,1:L])
    return m
end

function distill4(layer_list,L,q,p)
    CNOT = C"XX _X Z_ ZZ"
    # base case
    if q == 0
        return MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    end

    m1 = distill4(layer_list,L,q-1,p)
    m1 = gen_stabilizers(m1,layer_list[q],L)
    m1 = erasure_errors(m1,p)

    m2 = distill4(layer_list,L,q-1,p)
    m2 = gen_stabilizers(m2,layer_list[q],L)
    m2 = erasure_errors(m2,p)

    m3 = distill4(layer_list,L,q-1,p)
    m3 = gen_stabilizers(m3,layer_list[q],L)
    m3 = erasure_errors(m3,p)

    m4 = distill4(layer_list,L,q-1,p)
    m4 = gen_stabilizers(m4,layer_list[q],L)
    m4 = erasure_errors(m4,p)

    m5 = distill4(layer_list,L,q-1,p)
    m5 = gen_stabilizers(m5,layer_list[q],L)
    m5 = erasure_errors(m5,p)

    m = m1 ⊗ m2 ⊗ m3
    #println(rank(m))

    # check bit flip errors
    for i in 1:L
        apply!(m,CNOT,[i,i+L])
    end
    for i in 1:L
        project!(m,single_z(3*L,L+i),phases=false)
    end
    traceout!(m,(L+1):(2*L))

    #check phase flip errors
    for i in 1:L
        apply!(m,CNOT,[2*L+i,i])
    end
    for i in 1:L
        project!(m,single_x(3*L,2*L+i),phases=false)
    end
    traceout!(m,(2*L+1):(3*L))

    m = MixedDestabilizer(stabilizerview(m)[:,1:L])
    m = m ⊗ m4 ⊗ m5

    # check bit flip errors
    for i in 1:L
        apply!(m,CNOT,[i,i+L])
    end
    for i in 1:L
        project!(m,single_z(3*L,L+i),phases=false)
    end
    # traceout!(m,(L+1):(2*L))

    #check phase flip errors
    for i in 1:L
        apply!(m,CNOT,[2*L+i,i])
    end
    for i in 1:L
        project!(m,single_x(3*L,2*L+i),phases=false)
    end
    # traceout!(m,(2*L+1):(3*L))
    m = MixedDestabilizer(stabilizerview(m)[:,1:L])

    return m
end

function distill5(layer_list,L,q,p)
    CNOT = C"XX _X Z_ ZZ"
    # base case
    if q == 0
        input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
        m = gen_stabilizers(input,layer_list)
        println("m rank: $(rank(m)): $(stabilizerview(m))")
        m = erasure_errors(m,p)
        return m
    end

    m1 = distill5(layer_list,L,q-1,p)
    println("m1 after erasures rank: $(rank(m1)): $(stabilizerview(m1))")
    m2 = distill5(layer_list,L,q-1,p)
    println("m2 after erasures rank: $(rank(m2)): $(stabilizerview(m2))")

    m3 = distill5(layer_list,L,q-1,p)
    println("m3 after erasures rank: $(rank(m3)): $(stabilizerview(m3))")

    m = m1 ⊗ m2 ⊗ m3
    println("m concat rank: $(rank(m)): $(stabilizerview(m))")
    #println(rank(m))
    # check bit flip errors
    for i in 1:L
        apply!(m,CNOT,[i,i+L])
    end

    for i in 1:L
        project!(m,single_z(3*L,L+i),phases=false)
    end
    # traceout!(m,(L+1):(2*L))

    #check phase flip errors
    for i in 1:L
        apply!(m,CNOT,[2*L+i,i])
    end
    for i in 1:L
        project!(m,single_x(3*L,2*L+i),phases=false)
    end
    println("m after steane rank: $(rank(m)): $m")
    # traceout!(m,(2*L+1):(3*L))
    m = MixedDestabilizer(stabilizerview(m)[:,1:L])
    println("output m rank: $(rank(m)): $(stabilizerview(m))")
    return m
end

function channel_coding_cliff_east_non_rand()
    L_list = [51,52,53,54,55,56,57,58,59,60]
    p = 0.01
    num_samples = 10000
    open("data/cliff_east_non_rand_channel_coding.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for L in L_list
            for d in 1:20
                entropy_list = []
                for _ in 1:num_samples
                    #layer_list = gen_cliff_east_code_layers(L, d)
                    #randomize inputs
                    # xstab_ind = sample(1:L,L÷3,replace=false)
                    # zstab_ind = sample([i for i in 1:L if i ∉ xstab_ind],L÷3,replace=false)
                    # stabilizer_ind = vcat(xstab_ind,zstab_ind)
                    # input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))
                    # code = gen_stabilizers(input,layer_list)
                    code,stabilizer_ind = gen_cliff_east_code(L, d)
                    #println(code)
                    m = erasure_errors(copy(code),p)
                    # println("M")
                    # println(m)
                    measure_stabilizers(code,m,stabilizer_ind)
                    # println("after measurement")
                    # println(m)
                    entropy = (L - rank(m))/L
                    push!(entropy_list, entropy)
                end
                avg_entropy = sum(entropy_list)/num_samples
                std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
                sem_entropy = std_entropy/sqrt(num_samples)
                println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
                write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
            end
        end
    end
end
function channel_coding()
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/cliff_east_rand_input3_channel_coding.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 1:20
            entropy_list = []
            for _ in 1:num_samples
                layer_list = gen_cliff_east_code_layers(L, d)
                #randomize inputs
                xstab_ind = sample(1:L,L÷3,replace=false)
                zstab_ind = sample([i for i in 1:L if i ∉ xstab_ind],L÷3,replace=false)
                stabilizer_ind = vcat(xstab_ind,zstab_ind)
                input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))
                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding5()
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/cliff_east_rand_input5_channel_coding.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 1:20
            entropy_list = []
            for _ in 1:num_samples
                layer_list = gen_cliff_east_code_layers(L, d)
                #randomize inputs
                stabilizer_ind = sample(1:L,2*(L÷3),replace=false)
                sort!(stabilizer_ind)
                xstab_ind = [v for (i,v) in enumerate(stabilizer_ind) if i%2==1]
                input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))
                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding4()
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/cliff_east_rand_input4_channel_coding.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 1:20
            entropy_list = []
            for _ in 1:num_samples
                layer_list = gen_cliff_east_code_layers(L, d)
                #randomize inputs
                xstab_ind = []
                zstab_ind = []
                for i in 1:3:L
                    xi = rand(i:(i+2))
                    push!(xstab_ind,xi)
                    zi = rand([j for j in i:(i+2) if j != xi])
                    push!(zstab_ind,zi)
                end
                stabilizer_ind = vcat(xstab_ind,zstab_ind)
                input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))
                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding_block_wrapper()
    channel_coding11("cliffeast")
    channel_coding11("randcss1")
    channel_coding11("randcss2")
end
function channel_coding11(c)
    p = 0.01
    num_blocks = 8
    block_size=7
    L=num_blocks*block_size
    num_samples = 10000
    open("data/c_$(c)_channel_coding_block.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 1:10
            entropy_list = []
            for _ in 1:num_samples
                if c == "cliffeast"
                    layer_list = gen_cliff_east_code_open_block_layers(num_blocks,d,block_size)
                elseif c == "randcss1"
                    layer_list = gen_randcss_code_open_block_layers(num_blocks,d,block_size)
                elseif c == "randcss2"
                    layer_list = gen_randcss2_code_open_block_layers(num_blocks,d,block_size,0.5)
                end
                stabilizers = []
                midpoint = (block_size ÷ 2)+1
                shift = 0
                stabilizer_ind = []
                for _ in 1:num_blocks
                    for i in 1:(midpoint-1)
                        push!(stabilizer_ind,i+shift)
                        if i % 2 == 1
                            push!(stabilizers,single_x(L,i+shift))
                        else
                            push!(stabilizers,single_z(L,i+shift))
                        end
                    end
                    push!(stabilizers, single_z(L,midpoint+shift))
                    for i in (midpoint+1):block_size
                        push!(stabilizer_ind,i+shift)
                        if i % 2 == 1
                            push!(stabilizers,single_z(L,i+shift))
                        else
                            push!(stabilizers,single_x(L,i+shift))
                        end
                    end
                    shift+=block_size
                end
                input = MixedDestabilizer(Stabilizer([s for s in stabilizers]))
                #randomize inputs
                # xstab_ind = []
                # zstab_ind = []
                # for i in 1:3:L
                #     xi = rand(i:(i+2))
                #     push!(xstab_ind,xi)
                #     zi = rand([j for j in i:(i+2) if j != xi])
                #     push!(zstab_ind,zi)
                # end
                # stabilizer_ind = vcat(xstab_ind,zstab_ind)
                # input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))

                # input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
                # stabilizer_ind = [i for i in 1:L if (i+1)%3 != 0]

                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding6_wrapper()
    ph_list = [0.5,0.4,0.3,0.2,0.1,0.05]
    ph_str_list = ["0_5","0_4","0_3","0_2","0_1","0_05"]
    for (ph,ph_str) in zip(ph_list,ph_str_list)
        channel_coding6(ph,ph_str)
    end
end
function channel_coding6(ph,ph_str)
    println("ph: $(ph)")
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/rand_css_v5_ph_$(ph_str)_hxn_channel_coding_non_random_input.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 1:20
            entropy_list = []
            for _ in 1:num_samples
                layer_list = gen_css_code_v7_layers(L, d, ph)
                #randomize inputs
                # xstab_ind = []
                # zstab_ind = []
                # for i in 1:3:L
                #     xi = rand(i:(i+2))
                #     push!(xstab_ind,xi)
                #     zi = rand([j for j in i:(i+2) if j != xi])
                #     push!(zstab_ind,zi)
                # end
                # stabilizer_ind = vcat(xstab_ind,zstab_ind)
                # input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))

                input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
                stabilizer_ind = [i for i in 1:L if (i+1)%3 != 0]

                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding10(ph,ph_str)
    println("ph: $(ph)")
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/rand_css_v10_ph_$(ph_str)_channel_coding_non_random_input.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 1:20
            entropy_list = []
            for _ in 1:num_samples
                layer_list = gen_css_code_v10_layers(L, d, ph)
                #randomize inputs
                # xstab_ind = []
                # zstab_ind = []
                # for i in 1:3:L
                #     xi = rand(i:(i+2))
                #     push!(xstab_ind,xi)
                #     zi = rand([j for j in i:(i+2) if j != xi])
                #     push!(zstab_ind,zi)
                # end
                # stabilizer_ind = vcat(xstab_ind,zstab_ind)
                # input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))

                input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
                stabilizer_ind = [i for i in 1:L if (i+1)%3 != 0]

                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding9()
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/cliff_east_channel_coding_rand_input_p_0_01_L_50.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 1:20
            entropy_list = []
            for _ in 1:num_samples
                #layer_list = gen_css_code_v7_layers(L, d, 0.01)
                layer_list = gen_cliff_east_code_layers(L, d)
                #randomize inputs
                xstab_ind = []
                zstab_ind = []
                for i in 1:3:L
                    xi = rand(i:(i+2))
                    push!(xstab_ind,xi)
                    zi = rand([j for j in i:(i+2) if j != xi])
                    push!(zstab_ind,zi)
                end
                stabilizer_ind = vcat(xstab_ind,zstab_ind)
                input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))
                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding7(ph,ph_str)
    println("ph: $(ph)")
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/rand_css_v5_ph_$(ph_str)_channel_coding_mid_circuit_errors.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 11:15
            entropy_list = []
            for _ in 1:num_samples
                layer_list = gen_css_code_v5_layers(L, d, ph)
                #randomize inputs
                xstab_ind = []
                zstab_ind = []
                for i in 1:3:L
                    xi = rand(i:(i+2))
                    push!(xstab_ind,xi)
                    zi = rand([j for j in i:(i+2) if j != xi])
                    push!(zstab_ind,zi)
                end
                stabilizer_ind = vcat(xstab_ind,zstab_ind)
                input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))
                code = gen_stabilizers(copy(input),layer_list)
                m = mid_circuit_errors(copy(input),layer_list,p)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                #m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding8()
    L = 50
    p = 0.01
    num_samples = 10000
    open("data/test_cliff_east_channel_coding_mid_circuit_errors_p_0_01.csv", "a") do file
        write(file, "L,p,d,avg_entropy,sem_entropy\n")
        for d in 11:15
            entropy_list = []
            layer_list = gen_cliff_east_code_layers(L, d)
            for _ in 1:num_samples
                #randomize inputs
                xstab_ind = []
                zstab_ind = []
                for i in 1:3:L
                    xi = rand(i:(i+2))
                    push!(xstab_ind,xi)
                    zi = rand([j for j in i:(i+2) if j != xi])
                    push!(zstab_ind,zi)
                end
                stabilizer_ind = vcat(xstab_ind,zstab_ind)
                input = MixedDestabilizer(Stabilizer([ifelse(i in xstab_ind,single_x(L,i),single_z(L,i)) for i in 1:L]))
                code = gen_stabilizers(copy(input),layer_list)
                m = mid_circuit_errors(copy(input),layer_list,p)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                #m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function channel_coding3()
    L=48
    p = 0.01
    b=6
    num_samples = 10000
    open("data/cliff_east_block_L_$(L)_b_$(b)_channel_coding.csv", "a") do file
        write(file, "L,p,t,b,avg_entropy,sem_entropy\n")
        for d in 1:15
            entropy_list = []
            for _ in 1:num_samples
                layer_list = gen_cliff_east_code_periodic_block_layers(L, d, b)
                #randomize inputs
                input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
                stabilizer_ind = [i for i in 1:L if (i+1)%3 != 0]

                code = gen_stabilizers(input,layer_list)
                #code,stabilizer_ind = gen_cliff_east_code(L, d)
                #println(code)
                m = erasure_errors(copy(code),p)
                # println("M")
                # println(m)
                measure_stabilizers(code,m,stabilizer_ind)
                # println("after measurement")
                # println(m)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
            write(file, "$L,$p,$d,$b,$(avg_entropy),$(sem_entropy)\n")
        end
    end
end
function plot_mid_circuit_t_5()
    plot()
    filename1 = "data/mid_circuit_t_5_no_steane_ns_10000_L_50_p_0_001.csv"
    df1 = CSV.read(filename1, DataFrame)
    filename2 = "data/t_5_L_50_p_0_001_ns_50000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df1[!, "d"], df1[!, "avg_entropy"],yerr=df1[!, "sem_entropy"],label="no steane",title="mid circuit errors, L=50,p=0.001")
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="mid circuit errors, L=50,p=0.001")

    savefig("plots/mid_circuit_t_5_L_50_p_0_001.pdf")
end
function plot_mid_circuit_t_5_steane()
    plot()
    filename2 = "data/t_5_L_50_p_0_001_ns_50000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="mid circuit errors, L=50,p=0.001")

    savefig("plots/mid_circuit_w_steane_t_5_L_50_p_0_001.pdf")
end
function plot_t_layer_errors_t_5()
    plot()
    filename1 = "t_5_L_50_p_0_02_ns_5000_no_steane_errors_after_t_layers.csv"
    df1 = CSV.read(filename1, DataFrame)
    filename2 = "data/t_layer_errors_t_5_L_50_p_0_02_ns_50000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df1[!, "d"], df1[!, "avg_entropy"],yerr=df1[!, "sem_entropy"],label="no steane",title="errors every 5 layers, L=50,p=0.02")
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="errors every 5 layers, L=50,p=0.02")

    savefig("plots/t_layer_errors_t_5_L_50_p_0_02.pdf")
end
function plot_t_layer_errors_t_5_steane()
    plot()
    filename2 = "data/t_layer_errors_t_5_L_50_p_0_02_ns_50000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="errors every 5 layers, L=50,p=0.02")

    savefig("plots/t_layer_errors_w_steane_t_5_L_50_p_0_02.pdf")
end
function plot_mid_circuit_t_5_steane_v2()
    plot()
    filename2 = "data/t_5_L_50_p_0_01_ns_5000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="mid circuit errors, L=50,p=0.01")

    savefig("plots/mid_circuit_w_steane_t_5_L_50_p_0_01_ns_100_000_vmerged.pdf")
end
function plot_mid_circuit_t_5_steane_v3()
    plot()
    filename2 = "data/t_5_L_50_p_0_002_ns_5000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="mid circuit errors, L=50,p=0.002")

    savefig("plots/mid_circuit_w_steane_t_5_L_50_p_0_002_ns_100_000_vmerged.pdf")
end
function plot_mid_circuit_t_5_steane_v4()
    plot()
    filename2 = "data/t_5_L_50_p_0_003_ns_5000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="mid circuit errors, L=50,p=0.003")

    savefig("plots/mid_circuit_w_steane_t_5_L_50_p_0_003_ns_100_000_vmerged.pdf")
end
function plot_mid_circuit_t_5_steane_v5()
    plot()
    filename2 = "data/t_5_L_50_p_0_004_ns_5000_vmerged_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!(df2[!, "d"], df2[!, "avg_entropy"],yerr=df2[!, "sem_entropy"],label="with steane",title="mid circuit errors, L=50,p=0.004")

    savefig("plots/mid_circuit_w_steane_t_5_L_50_p_0_004_ns_100_000_vmerged.pdf")
end
function plot_mid_circuit_v4()
    plot()
    filename2 = "data/mid_circuit_errors_v4_c_cliffeast_t_5_L_50_p_0_001_ns_1000_vmerged_post_select2_n_5_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!([5,10,15], [df2[2, "avg_entropy"],df2[4, "avg_entropy"],df2[6, "avg_entropy"]],yerr=[df2[2, "sem_entropy"],df2[4, "sem_entropy"],df2[6, "sem_entropy"]],title="mid circuit errors, L=50,p=0.001,post-select=5")

    savefig("plots/mid_circuit_errors_v4_c_cliffeast_t_5_L_50_p_0_001_ns_1000_vmerged_post_select2_n_5_avg.pdf")
end
function plot_mid_circuit_v4_p_0_01()
    plot()
    filename2 = "data/mid_circuit_errors_v4_c_cliffeast_t_5_L_50_p_0_01_ns_1000_vmerged_post_select2_n_5_avg.csv"
    df2 = CSV.read(filename2, DataFrame)
    plot!([5,10,15], [df2[1, "avg_entropy"],df2[3, "avg_entropy"],df2[5, "avg_entropy"]],yerr=[df2[1, "sem_entropy"],df2[3, "sem_entropy"],df2[5, "sem_entropy"]],title="mid circuit errors, L=50,p=0.01,post-select=5")

    savefig("plots/mid_circuit_errors_v1_c_cliffeast_t_5_L_50_p_0_01_ns_1000_vmerged_post_select2_n_5_avg.pdf")
end
function plot_cnot_h()
    plot()
    #p_list = [0.5,0.4,0.3,0.2,0.1,0.08,0.06,0.04,0.02]
    #p_str_list = ["0_5","0_4","0_3","0_2","0_1","0_08","0_06","0_04","0_02"]
    p_list = [0.5,0.1,0.02]
    p_str_list = ["0_5","0_1","0_02"]
    for (p,p_str) in zip(p_list,p_str_list)
        filename = "data/rand_css_v5_ph_$(p_str)_channel_coding.csv"
        df = CSV.read(filename, DataFrame)
        plot!(df[!, "d"], df[!, "avg_entropy"],label="p=$p")
    end
    filename = "data/cliff_east_channel_coding_50_000.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],label="cliff-east")
    savefig("plots/channel_coding_cliff_east_plus_H_w_vary_p.pdf")
end
function plot_block_channel_coding()
    plot()
    filename = "data/c_randcss1_channel_coding_block.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="unif random-css",title="block_size=7,num_blocks=8,L=56")

    filename = "data/c_randcss2_channel_coding_block.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="CNOT+H random-css")

    filename = "data/c_cliffeast_channel_coding_block.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="clifford east")

    savefig("plots/block_channel_coding_vary_css_code.pdf")
end
function plot_cliffeast_randinput_channel_coding()
    plot()
    filename = "data/cliff_east_channel_coding_50_000.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="cliffeast")

    filename = "data/cliff_east_rand_input5_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="rand loc,X Z every other")

    filename = "data/cliff_east_rand_input4_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="rand loc, X Z log every 3")

    filename = "data/cliff_east_rand_input3_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="rand loc, X Z rand")

    savefig("plots/channel_coding_cliffeast_randinput.pdf")
end
function plot_cliffeast_defect_channel_coding()
    plot()
    filename = "data/cliff_east_channel_coding_50_000.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="cliffeast")

    filename = "data/rand_css_v10_ph_0_1_channel_coding_non_random_input.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="10% swaps at beginning")

    filename = "data/rand_css_v9_ph_0_1_channel_coding_non_random_input.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="10% unif rand css gates at beginning")

    filename = "data/rand_css_v8_ph_channel_coding_non_random_input.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="first layer replaced w/ unif rand css")

    filename = "data/rand_css_v7_ph_0_01_channel_coding_non_random_input.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="1% gates replaced w/ unif rand css")

    filename = "data/rand_css_v5_ph_0_1_hxn_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="10% layers followed by hadamards")
    filename = "data/rand_css_v5_ph_0_5_hxn_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="50% layers followed by hadamards")

    savefig("plots/channel_coding_cliffeast_defects.pdf")
end
function plot_retry_counts_randcss()
    plot()
    filename = "data/randcss_coding_block_retry_counts.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "t"], df[!, "avg_count"],yerr=df[!, "sem_count"])
    savefig("plots/retry_counts_randcss.pdf")
end
function plot_mid_circuit_errors_v5_steane_errors()
    plot()
    filename = "data/mid_circuit_errors_v5_c_cliffeast_ph_0_1_t_5_L_50_p_0_01_ns_5000_vmerged_post_select2_n_10_avg.csv"
    df = CSV.read(filename, DataFrame)
    plot!([df[2, "d"],df[4, "d"],df[6, "d"]], [df[2, "avg_entropy"],df[4, "avg_entropy"],df[6, "avg_entropy"]],yerr=[df[2, "sem_entropy"],df[4, "sem_entropy"],df[6, "sem_entropy"]],title="1% p,post-select=10")
    savefig("plots/mid_circuit_errors_v5_c_cliffeast_post_select_n_10_p_0_01.pdf")
end

function plot_channel_coding()
    plot()
    filename = "data/random_css_code_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="random-css")

    filename = "data/random_css_code_v2_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="random-css-CNOT+H")

    filename = "data/cliff_east_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="clifford-east")

    filename = "data/random_noncss_code_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="noncss",title="L=50,p=0.01",ylabel="avg entropy",xlabel="d")
    savefig("plots/channel_coding_vary_css_code.pdf")
end
function plot_channel_coding_v2()
    plot()
    filename = "data/cliff_east_channel_coding_50_000.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="clifford-east")

    filename = "data/rand_css_channel_coding_v3.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="unif rand css")

    filename = "data/rand_css2_channel_coding.csv"
    df = CSV.read(filename, DataFrame)
    plot!(df[!, "d"], df[!, "avg_entropy"],yerr=df[!, "sem_entropy"],label="rand CNOT+H")

    savefig("plots/channel_coding_vary_css_ns_10_000.pdf")
end
function plot_steane()
    round0 =  [0.007379999999999954,0.007039999999999952,0.004899999999999988,0.004439999999999995,0.002940000000000002,0.0019000000000000013,0.0008800000000000005,0.0012200000000000006,0.0011400000000000006,0.0008600000000000004]

    round1 = [0.00473999999999999,0.002700000000000002,0.0019400000000000012,0.0011800000000000007,0.00026,8.0e-5,0.00014000000000000001,5.9999999999999995e-5,2.0e-5,2.0e-5]
    round2 = [0.004719999999999988,0.002360000000000002,0.0012400000000000007,0.0011600000000000007,2.0e-5,4.0e-5,2.0e-5,0.0,0.0,0.0]

    round3 = [0.004339999999999997,0.0021600000000000013,0.0014000000000000009,0.0013600000000000007,8.0e-5,5.9999999999999995e-5,2.0e-5,0.0,5.9999999999999995e-5,0.0]
    plot()
    d = [i for i in 1:10]
    plot!(d,round0,label="no Stn Ver",title="L=50,p=0.02",xlabel="d",ylabel="avg entropy")
    plot!(d,round1,label="1 round of Stn Ver")
    plot!(d,round2,label="2 rounds of Stn Ver")
    plot!(d,round3,label="3 rounds of Stn Ver")
    savefig("plots/steane_ver_vary_rounds_no_mid_circuit_errors.pdf")
end
function plot_steane2()
    round1 = [0.00473999999999999,0.002700000000000002,0.0019400000000000012,0.0011800000000000007,0.00026,8.0e-5,0.00014000000000000001,5.9999999999999995e-5,2.0e-5,2.0e-5]
    round2 = [0.004719999999999988,0.002360000000000002,0.0012400000000000007,0.0011600000000000007,2.0e-5,4.0e-5,2.0e-5,0.0,0.0,0.0]

    round3 = [0.004339999999999997,0.0021600000000000013,0.0014000000000000009,0.0013600000000000007,8.0e-5,5.9999999999999995e-5,2.0e-5,0.0,5.9999999999999995e-5,0.0]
    plot()
    d = [i for i in 1:10]
    plot!(d[5:10],round1[5:10],label="1 round of Stn Ver",title="L=50,p=0.02",xlabel="d",ylabel="avg entropy")
    plot!(d[5:10],round2[5:10],label="2 rounds of Stn Ver")
    plot!(d[5:10],round3[5:10],label="3 rounds of Stn Ver")
    savefig("plots/steane_ver_vary_rounds_no_mid_circuit_errors_v2.pdf")
end
function channel_coding2()
    L = 50
    p = 0.02
    num_samples = 1000
    for d in 1:10
        entropy_list = []
        for _ in 1:num_samples
            layer_list = gen_cliff_east_code_layers(L, d)
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(input,layer_list,L)
            #code,stabilizer_ind = gen_cliff_east_code(L, d)
            #println(code)
            m = erasure_errors(copy(code),p)
            # println("M")
            # println(m)
            measure_stabilizers(code,m,stabilizer_ind)
            # println("after measurement")
            # println(m)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        println("d=$d avg_entropy: $(avg_entropy)")
    end
end
# function channel_coding()
#     L = 50
#     p = 0.05
#     num_samples = 1000
#     for d in 1:10
#         entropy_list = []
#         for _ in 1:num_samples
#             layer_list = gen_cliff_east_code_layers(L, d)
#             m = Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L])
#             stabilizer_ind = [i for i in 1:L if i%3 != 0]
#             code = gen_stabilizers(copy(m),layer_list,L)
#             #println(code)
#             m = erasure_errors(copy(code),p)
#             # println("M")
#             # println(m)
#             measure_stabilizers(code,m,stabilizer_ind)
#             # println("after measurement")
#             # println(m)
#             entropy = (L - rank(MixedStabilizer(m)))/L
#             push!(entropy_list, entropy)
#         end
#         avg_entropy = sum(entropy_list)/num_samples
#         println("avg_entropy: $(avg_entropy)")
#     end
# end
function mid_circuit_test()
    #m = Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L])
    #stabilizer_ind = [i for i in 1:L if i%3 != 0]

    L = 50
    p = 0.05
    num_samples = 1000
    for d in 1:10
        entropy_list = []
        for _ in 1:num_samples
            layer_list = gen_cliff_east_code_layers(L, d)
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(copy(input),layer_list,L)
            m = mid_circuit_errors(copy(input),layer_list,L,p)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        println("avg_entropy: $(avg_entropy)")
    end
end
function mid_circuit_distill_test()
    L = 20
    p = 0.001
    num_samples = 10000
    for d in 1:10
        entropy_list = []
        for _ in 1:num_samples
            layer_list = gen_cliff_east_code_layers(L, d)
            m = distill(layer_list,L,length(layer_list),p)
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(input,layer_list,L)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        println("avg_entropy: $(avg_entropy)")
    end
end

function mid_circuit_distill_test2()
    L = 50
    p = 0.02
    num_samples = 1000
    t=2
    for d in t:t:(3*t)
    #for d in [20]
        entropy_list = []
        for _ in 1:num_samples
            #layer_list = gen_cliff_east_code_layers(L, d)
            layer_list = gen_css_code_layers(L, d)
            i=1
            new_layer_list = []
            while i+t-1<=d
                push!(new_layer_list,layer_list[i:(i+t-1)])
                i+=t
            end
            m = distill2(new_layer_list,L,length(new_layer_list),p)
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(input,layer_list,L)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        println("avg_entropy: $(avg_entropy)")
    end
end
function mid_circuit_test2()
    L = 50
    p = 0.02
    num_samples = 5000
    t=5
    for d in t:t:(3*t)
    #for d in [20]
        entropy_list = []
        for _ in 1:num_samples
            #layer_list = gen_cliff_east_code_layers(L, d)
            layer_list = gen_cliff_east_code_layers(L, d)
            i=1
            new_layer_list = []
            while i+t-1<=d
                push!(new_layer_list,layer_list[i:(i+t-1)])
                i+=t
            end
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(copy(input),layer_list,L)
            m = copy(input)
            for layers in new_layer_list
                m = gen_stabilizers(m,layers,L)
                m = erasure_errors(m,p)
            end
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
        sem_entropy = std_entropy/sqrt(num_samples)
        println("avg_entropy: $(avg_entropy),sem_entropy: $(sem_entropy)")
    end
end

function mid_circuit_distill_test3()
    L = 50
    p = 0.0001
    num_samples = 50000
    t=5
    df = DataFrame(trial = 1:num_samples)
    for d in [5,10,15]
        println("d = $d")
        entropy_list = []
        layer_list = gen_cliff_east_code_layers(L, d)
        i=1
        new_layer_list = []
        while i+t-1<=d
            push!(new_layer_list,layer_list[i:(i+t-1)])
            i+=t
        end
        input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
        stabilizer_ind = [i for i in 1:L if i%3 != 0]
        code = gen_stabilizers(input,layer_list,L)
        for j in 1:num_samples
            m = distill3(new_layer_list,L,length(new_layer_list),p)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
            if j % 1000 == 0
                print("$(j)/$(num_samples): ")
                avg_entropy = sum(entropy_list)/j
                std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(j - 1))
                sem_entropy = std_entropy/sqrt(j)
                println("avg_entropy: $(avg_entropy),sem_entropy: $(sem_entropy)")
            end
        end
        df[!,"depth-$d"] = entropy_list
        println("DONE")
        avg_entropy = sum(entropy_list)/num_samples
        std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
        sem_entropy = std_entropy/sqrt(num_samples)
        println("avg_entropy: $(avg_entropy),sem_entropy: $(sem_entropy)")
    end
    CSV.write(filename,df)
end

function mid_circuit_test3()
    #m = Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L])
    #stabilizer_ind = [i for i in 1:L if i%3 != 0]

    L = 50
    p = 0.001
    num_samples = 10000
    t=5
    for d in t:t:(3*t)
        entropy_list = []
        for _ in 1:num_samples
            layer_list = gen_cliff_east_code_layers(L, d)
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(copy(input),layer_list,L)
            m = mid_circuit_errors(copy(input),layer_list,L,p)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
        sem_entropy = std_entropy/sqrt(num_samples)
        println("avg_entropy: $(avg_entropy),sem_entropy: $(sem_entropy)")
    end
end

function mid_circuit_distill_test4()
    L = 50
    p = 0.02
    num_samples = 1000
    t=5
    for d in t:t:(3*t)
    #for d in [20]
        entropy_list = []
        for _ in 1:num_samples
            #layer_list = gen_cliff_east_code_layers(L, d)
            layer_list = gen_css_code_layers(L, d)
            i=1
            new_layer_list = []
            while i+t-1<=d
                push!(new_layer_list,layer_list[i:(i+t-1)])
                i+=t
            end
            m = distill4(new_layer_list,L,length(new_layer_list),p)
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(input,layer_list,L)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        println("avg_entropy: $(avg_entropy)")
    end
end

function mid_circuit_distill_test5()
    L = 50
    p = 0.02
    num_samples = 1000
    for d in 1:10
    #for d in [20]
        entropy_list = []
        for _ in 1:num_samples
            layer_list = gen_cliff_east_code_layers(L, d)
            m = distill5(layer_list,L,3,p)
            input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
            stabilizer_ind = [i for i in 1:L if i%3 != 0]
            code = gen_stabilizers(input,layer_list,L)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        println("d=$d avg_entropy: $(avg_entropy)")
    end
end

function mid_circuit_distill_test5_v2()
    L = 50
    p = 0.01
    num_samples = 1
    for q in [1]
        println("q: $q")
        for d in [5]
            println("d: $d")
        #for d in [20]
            entropy_list = []
            for _ in 1:num_samples
                #layer_list = gen_cliff_east_code_layers(L, d)
                layer_list = gen_css_code_v5_layers(L, d, 0.1)
                m = distill5(layer_list,L,q,p)
                input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
                stabilizer_ind = [i for i in 1:L if i%3 != 0]
                code = gen_stabilizers(input,layer_list)
                measure_stabilizers(code,m,stabilizer_ind)
                entropy = (L - rank(m))/L
                push!(entropy_list, entropy)
            end
            avg_entropy = sum(entropy_list)/num_samples
            std_entropy = sqrt(sum((entropy_list .- avg_entropy).^2)/(num_samples - 1))
            sem_entropy = std_entropy/sqrt(num_samples)
            println("q=$q d=$d avg_entropy: $(avg_entropy), sem_entropy: $(sem_entropy)")
        end
    end
end


function mid_circuit_distill_test6()
    L = 50
    p = 0.001
    num_samples = 5000
    layer_list = gen_cliff_east_code_layers(L, 20)
    new_layer_list = [layer_list[1:5]]

    input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
    stabilizer_ind = [i for i in 1:L if i%3 != 0]
    code = gen_stabilizers(input,layer_list[1:5],L)
    entropy_list = []
    for _ in 1:num_samples
        m = distill3(new_layer_list,L,length(new_layer_list),p)
        measure_stabilizers(code,m,stabilizer_ind)
        entropy = (L - rank(m))/L
        push!(entropy_list, entropy)
    end
    avg_entropy = sum(entropy_list)/num_samples
    println("d=5 avg_entropy: $(avg_entropy)")

    for d in [6,7,8,9,10,11,12]
    #for d in [20]
        new_layer_list = [layer_list[1:5]]
        entropy_list = []
        push!(new_layer_list,layer_list[6:d])
        println("num rounds: $(length(new_layer_list))")
        input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
        stabilizer_ind = [i for i in 1:L if i%3 != 0]
        code = gen_stabilizers(input,layer_list[1:d],L)
        for _ in 1:num_samples
            m = distill3(new_layer_list,L,length(new_layer_list),p)
            measure_stabilizers(code,m,stabilizer_ind)
            entropy = (L - rank(m))/L
            push!(entropy_list, entropy)
        end
        avg_entropy = sum(entropy_list)/num_samples
        println("d=$d avg_entropy: $(avg_entropy)")
    end
end
#
# function mid_circuit_test2()
#     L = 50
#     p = 0.02
#     num_samples = 1000
#     t=5
#     for d in t:t:(3*t)
#     #for d in [20]
#         entropy_list = []
#         for _ in 1:num_samples
#             layer_list = gen_cliff_east_code_layers(L, d)
#             i=1
#             new_layer_list = []
#             while i+t-1<=d
#                 push!(new_layer_list,layer_list[i:(i+t-1)])
#                 i+=t
#             end
#             #m = distill2(new_layer_list,L,length(new_layer_list),p)
#             input = MixedDestabilizer(Stabilizer([ifelse(i%3==1,single_x(L,i),single_z(L,i)) for i in 1:L]))
#             stabilizer_ind = [i for i in 1:L if i%3 != 0]
#             code = gen_stabilizers(copy(input),layer_list,L)
#             m = copy(input)
#             for layers in new_layer_list
#                 m = gen_stabilizers(m,layers,L)
#                 m = erasure_errors(m,p)
#             end
#             measure_stabilizers(code,m,stabilizer_ind)
#             entropy = (L - rank(m))/L
#             push!(entropy_list, entropy)
#         end
#         avg_entropy = sum(entropy_list)/num_samples
#         println("avg_entropy: $(avg_entropy)")
#     end
# end
