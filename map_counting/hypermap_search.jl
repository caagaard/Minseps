include("perm_utils.jl")
include("combinadic.jl")
include("murnaghan_nakayam.jl")
using Oscar
using Combinatorics
#using FLoops
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

# Searches a chunk of specified conjugacy class for alphas that gives a hypermap of the desired genus
# where phi has no fixed points.
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
			# Check if (p_inv,alpha) is the canonical representative of its isomorphism class
                    for g in HH
                        if (alpha^g).d < alpha.d
                            is_min=0
                            break
                        end
                    end
                    if is_min ==1
			# If alpha is conjugate to p_inv we now check if (p_inv, alpha) is the canonical
			# representative of its `extended isomorphism class' which allows for changing the
			# 2-coloring of the vertices
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

# Searches the conjugacy class in S_n corresponding to part for values of alpha that gives a 
# hypermap of desired genus with phi having no fixed points
# Returns a list with a single representative of each `extended isomorphism class', where we
# extend the regular definition of hypermap isomorphism to include swapping the 2-coloring 
# of the vertices
function get_alphas_dist(n::Int, g::Int, sigmatemp::Vector{Vector{Int}}, part::Vector{Int}, n_theta_cycles::Int)
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
        sigma = cperm(S, sigmatemp...)
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
				# Check if (sigma, alpha, phi) is the canonical representative
				# of its isomorphism class
		            for g in HH
		                if (theta^g).d < theta.d
		                    is_min=0
		                    break
		                end
		            end
                    if is_min ==1
			# If sigma is conjugate to alpha, we check it it's the canonical representative
			# of the extended isomorphism class
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
function get_phi_candidates_v1(n::Int, part::Vector{Int}, g::Int, sigmatemp::Vector{Vector{Int}}, n_phi_cycles::Int, num_nodes::Int)
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
        threadload = div(total_load, (num_nodes*Threads.nthreads()))
		S = symmetric_group(n)
		sigma = cperm(S,sigmatemp...)
		p_inv = Perm(Vector{Int}(sigma^(-1)))
		H = centralizer(S,sigma)
		HH = [Perm(Vector{Int}(x)) for x in H[1]]
		outlist = [Vector{Int}[] for i in 1:Threads.nthreads()]
        for i in 1:Threads.nthreads()
            sizehint!(outlist[i], threadload)
        end
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

function get_phi_candidates_v1(n::Int, part::Vector{Int}, g::Int, sigmatemp::Vector{Vector{Int}}, n_phi_cycles::Int, num_nodes::Int, nm_flag::Int)
	if sum(part) != n
		println("Invalid Partition")
		return([])
	else
        cc= counter(part)
        K = reverse(sort([k for k in keys(cc)]))
        PP = perm_counter(vcat([[k for i in 1:cc[k]] for k in K]...))
        if n_phi_cycles ==1
            hintload = div(map_bound(n, part), n)
        elseif n_phi_cycles ==2
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
        flush(stdout)
	    S = symmetric_group(n)
	    sigma = cperm(S,sigmatemp...)
	    p_inv = Perm(Vector{Int}(sigma^(-1)))
	    H = centralizer(S,sigma)
	    HH = [Perm(Vector{Int}(x)) for x in H[1]]
	    outlist = [Vector{Int}[] for i in 1:Threads.nthreads()]
        for i in 1:Threads.nthreads()
            sizehint!(outlist[i], hintload)
        end
        Threads.@threads for i in 1:Threads.nthreads()
            if i<Threads.nthreads()
                find_phis(i, threadload, threadload, part, cc, K, p_inv, HH, PP, n_phi_cycles, n, outlist[i])
            else
                tload = total_load - threadload*(Threads.nthreads()-1)
                find_phis(i, threadload, tload, part, cc, K, p_inv, HH, PP, n_phi_cycles,n,outlist[i])
            end
        end
		return([[Vector{Int}(sigma), x] for x in reduce(vcat,outlist)])
	end
end


