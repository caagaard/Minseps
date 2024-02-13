include("perm_utils.jl")
include("combinadic.jl")
using Oscar
using Combinatorics
using FLoops
using DataStructures

function perm_components(part::Array{Array{Int,1},1}, n::Int)
        starters = [chunk[1] for chunk in part]
        perm_parts = [chunk[2:length(chunk)] for chunk in part]
        return(starters, perm_parts,n)
end
function perm_counter(s::Array{Int,1})
        return([1:factorial(i-1) for i in s])
end

function make_perm(starters::Array{Int,1}, perm_parts::Array{Array{Int,1},1}, n::Int, index::NTuple)
        phi = [Int(1):n;]
        for i in 1:length(starters)
                tempperm = vcat(starters[i], nthperm(perm_parts[i], index[i]))
                for j in 1:(length(tempperm)-1)
                        phi[tempperm[j]] = tempperm[j+1]
                        phi[tempperm[end]] = tempperm[1]
                end
        end
        return(Perm(phi))
end

function make_cartesian(v::Vector{Int})
       n = length(v) 
       return(CartesianIndices(ntuple(i-> 1:v[i], n)))
end

function get_phi_candidates_v2(n::Int,k::Int,g::Int, psitemp::Vector{Vector{Int}})
	S = symmetric_group(n)
	sigma = cperm(S, psitemp...)
	psi = Perm(Vector{Int}(sigma))
	needed_verts = n-k -length(cycles(psi))+2 - 2*g
	H = centralizer(S,sigma)
	HH = [Perm(Vector{Int}(x)) for x in H[1]]
        parts = [part for part in Combinatorics.partitions([Int(1):Int(n);],k)]
	outlist = [Perm{Int}[] for i in 1:Threads.nthreads()]
        for part in parts
                decomp = perm_components(part,Int(n))
                #PP = perm_counter([length(part[i]) for i in 1:length(part)])
                flooptime = time()
                #@floop for index in Iterators.product(PP...)
                PP = make_cartesian([factorial(length(part[i])-1) for i in 1:length(part)])
                Threads.@threads for index in PP
			phi = make_perm(decomp[1],decomp[2], decomp[3], Tuple.(index))
			theta = psi*phi
			if length(cycles(theta)) == needed_verts
				if 1 in [length(cyc) for cyc in cycles(theta)]
					nothing
				else
					is_min =1
					for g in HH
						if (g^(-1)*phi*g).d < phi.d
							is_min =0
							break
						end
					end
					if is_min ==1
						push!(outlist[Threads.threadid()], phi)
					end
				end
			end
                end     
        end     
        return([[Perm(Vector{Int}(sigma)), x] for x in reduce(vcat,outlist)])
end 

function find_phis(i::Int, avgtload::Int, tload::Int, part::Vector{Int}, cc::Accumulator{Int, Int}, K::Vector{Int}, p_inv::Perm{Int}, HH::Vector{Perm{Int}}, PP::Vector{UnitRange{Int}}, n_phi_cycles::Int, n::Int, outlist::Vector{Perm{Int}})
        # Need correct start position here, but this will allow for testing
        #inner_struct =  initialize(i)
        #iter_length = makelength(i)
        # need a "next" for conj class element
        #outlist = Perm{}[]
        combo_part = [[1:p;] for p in part]
        j=0
        while j<avgtload*(i-1)
                conj_class_next!(n,combo_part, [cc[k] for k in K])
                j = j+1
        end
        while j<tload+avgtload*(i-1)
                decomp = perm_components(ct_to_p(combo_part), n)
                for index in Iterators.product(PP...)
                        theta = make_perm(decomp[1], decomp[2], decomp[3], index)
                        phi = theta^(-1)*p_inv
                        if length(cycles(phi)) == n_phi_cycles
                                is_min =1

		                        for g in HH
		                                if (theta^g).d < theta.d
		                                        is_min=0
		                                break
		                                end
		                        end
                                if is_min ==1
                                        push!(outlist, phi)
                                end
                        end
                end
                combo_part, iter_done = conj_class_next!(n, combo_part, [cc[k] for k in K])
                j = j+1
        end
        #return(outlist)
end

function get_phi_candidates_v1(n::Int, part::Vector{Int}, g::Int, psitemp::Vector{Vector{Int}}, n_phi_cycles::Int)
	if sum(part) != n
		println("Invalid Partition")
		return([])
	else
        cc= counter(part)
        K = reverse(sort([k for k in keys(cc)]))
        PP = perm_counter(vcat([[k for i in 1:cc[k]] for k in K]...))
        total_load_denom = 1
        for k in K
                total_load_denom = total_load_denom*factorial(k)^(cc[k])*factorial(cc[k])
        end
        total_load = div(factorial(n), total_load_denom)
        threadload = div(total_load, Threads.nthreads())
		S = symmetric_group(n)
		sigma = cperm(S,psitemp...)
		p_inv = Perm(Vector{Int}(sigma^(-1)))
		H = centralizer(S,sigma)
		HH = [Perm(Vector{Int}(x)) for x in H[1]]
		outlist = [Perm{Int}[] for i in 1:Threads.nthreads()]
        for i in 1:Threads.nthreads()
                sizehint!(outlist[i], threadload)
        end
        #flooptime = time()
        Threads.@threads for i in 1:Threads.nthreads()
                if i<Threads.nthreads()
                        find_phis(i, threadload, threadload, part, cc, K, p_inv, HH, PP, n_phi_cycles, n, outlist[i])
                else
                        tload = total_load - threadload*(Threads.nthreads()-1)
                        find_phis(i, threadload, tload, part, cc, K, p_inv, HH, PP, n_phi_cycles,n,outlist[i])
                end
                #outlist[i] = Threads.@spawn find_phis(i,tload, part, cc, K, p_inv, HH, PP, n_phi_cycles, n)
        end
        #@sync begin
                #for i in 1:Threads.nthreads()
                        #outlist[1] = Threads.@spawn find_phis(i, tload, part, cc, K, p_inv, HH, PP, n_phi_cycles,n)
                #end
        #end
        #outlist[1] = find_phis(1, total_load, part,cc,K,p_inv, HH, PP, n_phi_cycles,n)
		return([[Perm(Vector{Int}(sigma)), x] for x in reduce(vcat,outlist)])
	end
end

function get_phi_candidates_thread(n::Int, part::Vector{Int}, g::Int, psitemp::Vector{Vector{Int}}, n_phi_cycles::Int)
	if sum(part) != n
		println("Invalid Partition")
		return([])
	else
		S = symmetric_group(n)
		sigma = cperm(S,psitemp...)
		p_inv = Perm(Vector{Int}(sigma^(-1)))
		H = centralizer(S,sigma)
		HH = [Perm(Vector{Int}(x)) for x in H[1]]
		outlist = [Perm{Int}[] for i in 1:Threads.nthreads()]
        cc= counter(part)
        blocklengths=[k*cc[k] for k in keys(cc)]
        K = [k for k in keys(cc)]
        outer_structure = [colex_bitstring(n - sum(blocklengths[1:i-1]), blocklengths[i]) for i in 1:length(blocklengths)]
        PP = perm_counter(vcat([[k for i in 1:cc[k]] for k in keys(cc)]...))
        looptime = time()
        #tuple_prod = Iterators.product([1:makeRCI(k,cc[k]).Length for k in keys(cc)]..., PP...)
        tuple_prod = Iterators.product(outer_structure..., PP...)
        v_coord = make_cartesian([makeRCI(k,cc[k]).Length for k in keys(cc)])
        #tuple_prod = Iterators.product(outer_structure..., [1:makeRCI(k,cc[k]).Length for k in keys(cc)]..., PP...)
        #Threads.@threads for vindex in collect(Iterators.product([1:makeRCI(k,cc[k]).Length for k in keys(cc)]...))
        Threads.@threads for vindex in v_coord
        #@floop for thetatuple in tuple_prod
            # thetatuple is a tuple with the first length(cc) being the big blocks and the remaining elements the "regular_combination_tuples"
            #block_parts = deepcopy(thetatuple[1:length(keys(cc))])
                v_parts = deepcopy([unrank_reg_combo(vindex[i], blocklengths[i], K[i]) for i in 1:length(vindex)])
            for thetatuple in tuple_prod
                #vindex = thetatuple[1:(end-length(part))]
                block_parts = thetatuple[1:length(keys(cc))]
                index = thetatuple[(1+end- length(part)):end]
                N = [1:n;]
                t_parts = Vector{Vector{Int}}[]
                for i in 1:length(keys(cc))
                        push!(t_parts, [[N[filter(x-> block_parts[i][x] == 1, eachindex(N))][j] for j in w] for w in ct_to_p(v_parts[i])])
                        N = N[filter(x -> block_parts[i][x] == 0, eachindex(block_parts[i]))]
                end
                theta_part = vcat(t_parts...)
                decomp = deepcopy(perm_components(theta_part, n))
                theta = make_perm(decomp[1],decomp[2], decomp[3], index)
                phi = theta^(-1)*p_inv
                if length(cycles(phi)) == n_phi_cycles
                        is_min =1
		                for g in HH
		                        if (theta^g).d < theta.d
		                                is_min=0
		                                break
		                        end
		                end
		                if is_min ==1
		                        push!(outlist[Threads.threadid()], phi)
		                end
		        end
		    end
        end
		return([[Perm(Vector{Int}(sigma)), x] for x in reduce(vcat,outlist)])
	end
end


