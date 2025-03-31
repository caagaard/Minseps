include("perm_utils.jl")
include("combinadic.jl")
include("murnaghan_nakayam.jl")
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

# This needs to be made distributed somehow?
# When we ran it all cases that use v2 could run 
# easily on a single node
# Also need to change output for `local iso testing'

# Currently changing to keep old hypermaps output vs new graph output, just adjusting flow to 
# be more distributable.  To convert to graphs will need to add a new arraw of counters (one
# counter per thread)
function get_phi_candidates_v2(n::Int,k::Int,g::Int, psitemp::Vector{Vector{Int}})
    #println("running v2")
	#flush(stdout)
	S = symmetric_group(n)
	sigma = cperm(S, psitemp...)
	psi = Perm(Vector{Int}(sigma))
	needed_verts = n-k -length(cycles(psi))+2 - 2*g
	H = centralizer(S,sigma)
	HH = [Perm(Vector{Int}(x)) for x in H[1]]
    parts = [part for part in Combinatorics.partitions([Int(1):Int(n);],k)]
	outlist = [Vector{Int}[] for i in 1:Threads.nthreads()]
    # Need something more easily distributed here.  
    # Currently iterating over approximately all (reduced a bit from symmetry) permutations
    # with given number of cycles
    # Can either make unrank/next methods for this class or go over conj classes and use
    # conjclass methods
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
					    push!(outlist[Threads.threadid()], phi.d)
					end
				end
			end
        end     
    end     
    return([[Vector{Int}(sigma), x] for x in reduce(vcat,outlist)])
end 

function find_alphas(i::Int, avgtload::Int, tload::Int, par::Vector{Int}, cc::Accumulator{Int, Int}, K::Vector{Int}, p_inv::Perm{Int}, HH::Vector{Perm{Int}}, PP::Vector{UnitRange{Int}},
n_theta_cycles::Int, n::Int, outlist::Vector{Vector{Int}})
    combo_part = unrank_combo_partition(1+avgtload*(i-1), n, K, [cc[k] for k in K])
    for j in 1:tload
        decomp = perm_components(ct_to_p(combo_part), n)
        for index in Iterators.product(PP...)
            alpha = make_perm(decomp[1], decomp[2], decomp[3], index)
            theta = p_inv*alpha^(-1)
            if length(cycles(theta)) == n_theta_cycles
                if 1 in [length(cyc) for cyc in cycles(theta)]
                    nothing
                else
                    is_min =1
                    for g in HH
                        if (alpha^g).d < alpha.d
                            is_min=0
                            break
                        end
                    end
                    if is_min ==1
                        if conjclass(p_inv) == conjclass(alpha)
                            rho = solve_conjugation(alpha, p_inv^(-1),n)
                            pprime = rho^(-1)*(p_inv^(-1))*rho
                            for g in HH
                                if (pprime^g).d < alpha.d
                                    is_min =0
                                    break
                                end
                            end
                            if is_min ==1
                                push!(outlist, alpha.d)
                            end
                        else
                            push!(outlist, alpha.d)
                        end
                    end
                end
            end
        end
        combo_part, iter_done = conj_class_next!(n, combo_part, [cc[k] for k in K])
    end
end

function get_alphas_dist(n::Int, g::Int, psitemp::Vector{Vector{Int}}, part::Vector{Int}, n_theta_cycles::Int)
    #println("running alpha search")
    #flush(stdout)
    num_nodes =1
    if sum(part) != n
        println("Invalid Partition")
        return([])
    else
        cc = counter(part)
        K = reverse(sort([k for k in keys(cc)]))
        PP = perm_counter(vcat([[k for i in 1:cc[k]] for k in K]...))
        total_load_denom = 1
        for k in K
            total_load_denom = total_load_denom*factorial(k)^(cc[k])*factorial(cc[k])
        end
        total_load = div(factorial(n), total_load_denom)
        threadload = div(total_load, (num_nodes*Threads.nthreads()))
        S = symmetric_group(n)
        sigma = cperm(S, psitemp...)
        p_inv = Perm(Vector{Int}(sigma^(-1)))
        H = centralizer(S,sigma)
        HH=[Perm(Vector{Int}(x)) for x in H[1]]
        outlist = [Vector{Int}[] for i in 1:Threads.nthreads()]
        for i in 1:Threads.nthreads()
            sizehint!(outlist[i], threadload)
        end
        Threads.@threads for i in 1:Threads.nthreads()
            if i<Threads.nthreads()
                find_alphas(i, threadload, threadload, part, cc, K, p_inv, HH, PP, n_theta_cycles, n, outlist[i])
            else
                tload = total_load - threadload*(num_nodes*Threads.nthreads()-1)
                find_alphas(i, threadload, tload, part, cc, K, p_inv, HH, PP, n_theta_cycles,n,outlist[i])
            end
        end
        return([[Vector{Int}(sigma),x] for x in reduce(vcat,outlist)])
    end
end

# This function should run on each thread on each node
# Distributed version should not store hypermaps, just graphs
# This means ''adjusting count'' for color swaps
# Output should be list of graphs and the (color swap adjusted) count of hypermaps
function find_phis(i::Int, avgtload::Int, tload::Int, part::Vector{Int}, cc::Accumulator{Int, Int}, K::Vector{Int}, p_inv::Perm{Int}, HH::Vector{Perm{Int}}, PP::Vector{UnitRange{Int}}, n_phi_cycles::Int, n::Int, outlist::Vector{Vector{Int}})
    combo_part = unrank_combo_partition(1+avgtload*(i-1), n, K, [cc[k] for k in K])
    for j in 1:tload 
        decomp = perm_components(ct_to_p(combo_part), n)
        for index in Iterators.product(PP...)
            theta = make_perm(decomp[1], decomp[2], decomp[3], index)
            phi = theta^(-1)*p_inv
            if length(cycles(phi)) == n_phi_cycles
                if n_phi_cycles != length(cycles(p_inv)) || conjclass(phi)>= conjclass(p_inv)
                    is_min =1
		            for g in HH
		                if (theta^g).d < theta.d
		                    is_min=0
		                    break
		                end
		            end
                    if is_min ==1
                        if conjclass(phi) == conjclass(p_inv)
                            rho = solve_conjugation(phi, p_inv^(-1),n)
                            phip = rho^(-1)*(p_inv)*rho
                            thetap = p_inv*phip
                            for g in HH
                                if (thetap^g).d < theta.d
                                    is_min=0
                                    break
                                end
                            end
                            if is_min ==1
                                push!(outlist, phi.d)
                            end
                        else
                            push!(outlist, phi.d)
                        end
                    end
                end
            end
        end
        combo_part, iter_done = conj_class_next!(n, combo_part, [cc[k] for k in K])
    end
end

# This version will should run when memory bounding isn't needed
# Should return count of hypermaps (adjusted to account for color swap) and list of graphs
# Assumed that all nodes have same number of threads and same memory.
# If this is a problem we can eliminate multithreading and run fully distributed but
# this may cause performance loss for graph iso
# Should add a parameter for number of nodes this is num_nodes
function get_phi_candidates_v1(n::Int, part::Vector{Int}, g::Int, psitemp::Vector{Vector{Int}}, n_phi_cycles::Int, num_nodes::Int)
    #println("running v1 no nmb")
    #flush(stdout)
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
        # Should be adjusted to account for worker numbers
        total_load = div(factorial(n), total_load_denom)
        threadload = div(total_load, (num_nodes*Threads.nthreads()))
        #threadload = 18381
        #println(threadload)
        #flush(stdout)
		S = symmetric_group(n)
		sigma = cperm(S,psitemp...)
		p_inv = Perm(Vector{Int}(sigma^(-1)))
		H = centralizer(S,sigma)
		HH = [Perm(Vector{Int}(x)) for x in H[1]]
        # outlist needs to be updated to be a list of graphs
		outlist = [Vector{Int}[] for i in 1:Threads.nthreads()]
        for i in 1:Threads.nthreads()
            sizehint!(outlist[i], threadload)
        end
        #flooptime = time()
        #counters = zeros(Int64, Threads.nthreads())
        # find_phis needs worker id along with threadid to find start point
        # This loop should check worker id in the conditional (only last worker needs to use it)
        Threads.@threads for i in 1:Threads.nthreads()
            if i<Threads.nthreads()
                find_phis(i, threadload, threadload, part, cc, K, p_inv, HH, PP, n_phi_cycles, n, outlist[i])
            else
                tload = total_load - threadload*(num_nodes*Threads.nthreads()-1)
                find_phis(i, threadload, tload, part, cc, K, p_inv, HH, PP, n_phi_cycles,n,outlist[i])
            end
        end
		return([[Vector{Int}(sigma), x] for x in reduce(vcat,outlist)])
	end
end

function get_phi_candidates_v1(n::Int, part::Vector{Int}, g::Int, psitemp::Vector{Vector{Int}}, n_phi_cycles::Int, num_nodes::Int, nm_flag::Int)
	#println("running v1 with nmb")
	#flush(stdout)
	if sum(part) != n
		println("Invalid Partition")
		return([])
	else
        cc= counter(part)
        K = reverse(sort([k for k in keys(cc)]))
        PP = perm_counter(vcat([[k for i in 1:cc[k]] for k in K]...))
        if n_phi_cycles ==1
#       nmb = map_bound(n, part)
            hintload = div(map_bound(n, part), n)
        elseif n_phi_cycles ==2
#           nmb=0
            hintload = div(round(sum([map_bound(n, [n-i, i], part) for i in 1:div(n,2)])), n)
        else
            println("use version without mn rule")
            flush(stdout)
            return([])
        end
        total_load_denom = 1
        for k in K
            total_load_denom = total_load_denom*factorial(k)^(cc[k])*factorial(cc[k])
        end
        total_load = div(factorial(n), total_load_denom)
        threadload = div(total_load, Threads.nthreads())
        #println(threadload)
        flush(stdout)
	    S = symmetric_group(n)
	    sigma = cperm(S,psitemp...)
	    p_inv = Perm(Vector{Int}(sigma^(-1)))
	    H = centralizer(S,sigma)
	    HH = [Perm(Vector{Int}(x)) for x in H[1]]
	    outlist = [Vector{Int}[] for i in 1:Threads.nthreads()]
        #counters = zeros(Int64, Threads.nthreads())
        for i in 1:Threads.nthreads()
            sizehint!(outlist[i], hintload)
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
		return([[Vector{Int}(sigma), x] for x in reduce(vcat,outlist)])
	end
end


