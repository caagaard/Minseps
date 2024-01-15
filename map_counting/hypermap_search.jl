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

function get_phi_candidates_v2(n::Int,k::Int,g::Int, psitemp::Vector{Vector{Int}})
	S = symmetric_group(n)
	sigma = cperm(S, psitemp...)
	psi = Perm(Vector{Int}(sigma))
	needed_verts = n-k -length(cycles(psi))+2 - 2*g
	#println(psi)
	H = centralizer(S,sigma)
	HH = [Perm(Vector{Int}(x)) for x in H[1]]
        parts = [part for part in Combinatorics.partitions([Int(1):Int(n);],k)]
	outlist = []
	for i in 1:Threads.nthreads()
		push!(outlist, [])
	end
        for part in parts
                decomp = perm_components(part,Int(n))
                PP = perm_counter([length(x) for x in part])
                flooptime = time()
                @floop for index in Iterators.product(PP...)
                        phi = make_perm(decomp[1],decomp[2], decomp[3], index)
			            is_min =1
			            for g in HH
				                if (g^(-1)*phi*g).d < phi.d
					                    is_min =0
					                    break
				                end
			            end
			            if is_min ==1
				                theta = psi*phi
				                # Temporary version just for genus 3, n = 12, k=1
				                        if length(cycles(theta)) == needed_verts
					                            if 1 in [length(cyc) for cyc in cycles(theta)]
						                                nothing
					                            else
                        			                    push!(outlist[Threads.threadid()], phi)
					                            end
				                        end
			            end
                end     
                #println("floop time = ")
                #println(string(time()-flooptime))
        end     
        return([Perm(Vector{Int}(sigma)), reduce(vcat,outlist)])
end 

function get_phi_candidates_v1(n::Int, part::Vector{Int}, g::Int, psitemp::Vector{Vector{Int}}, n_phi_cycles::Int)
        #println(string(part))
        #flush(stdout)
	if sum(part) != n
		println("Invalid Partition")
		return([])
	else
		S = symmetric_group(n)
		sigma = cperm(S,psitemp...)
		psi = Perm(Vector{Int}(sigma))
		H = centralizer(S,sigma)
		HH = [Perm(Vector{Int}(x)) for x in H[1]]
		outlist = []
		for i in 1:Threads.nthreads()
			push!(outlist, [])
		end

        cc= counter(part)
        blocklengths=[k*cc[k] for k in keys(cc)]
        outer_structure = [colex_bitstring(n - sum(blocklengths[1:i-1]), blocklengths[i]) for i in 1:length(blocklengths)]
        #PP = perm_counter([x for x in sort(part)])
        #PP = perm_counter([x for x in part])
        PP = perm_counter(vcat([[k for i in 1:cc[k]] for k in keys(cc)]...))
		#thetap = Perm(Vector{Int}(cperm(S, make_default_perm(part)...)))
        #thetap = cperm(S, make_default_perm(part)...)
		#CC = [Perm(Vector{Int}(x)) for x in conjugacy_class(S, thetap)]
        #CC = 
        flooptime = time()
		tuple_prod= Iterators.product(outer_structure..., [makeRCI(k, cc[k]) for k in keys(cc)]..., PP...)
        #tuple_prod= Iterators.product(outer_structure..., [makeRCI(k, cc[k]) for k in keys(cc)]...)
        @floop for thetatuple in tuple_prod
            # thetatuple is a tuple with the first length(cc) being the big blocks and the remaining elements the "regular_combination_tuples"
            # we need 
            block_parts = deepcopy(thetatuple[1:length(keys(cc))])
            #println(block_parts)
            v_parts = deepcopy([x for x in thetatuple[(length(keys(cc))+1):(end-length(part))]])
            #v_parts = deepcopy([x for x in thetatuple[(length(keys(cc))+1):end]])
            index = thetatuple[(1+end- length(part)):end]
            #println(index)
            #println(v_parts)
            N = [1:n;]
            t_parts = Vector{Vector{Int}}[]
            for i in 1:length(keys(cc))
                push!(t_parts, [[N[filter(x-> block_parts[i][x] == 1, eachindex(N))][j] for j in w] for w in ct_to_p(v_parts[i])])
                N = N[filter(x -> block_parts[i][x] == 0, eachindex(block_parts[i]))]
            end
            theta_part = vcat(t_parts...)
            #print("checking lengths match: ")
            #println(string(length(theta_part == length(part))))
            #println(string(theta_part))
            #println(string(part))
            #println(index)
            decomp = deepcopy(perm_components(theta_part, n))
            #PP = perm_counter([length(theta_part[i]) for i in 1:length(theta_part)])
            #for index in Iterators.product(PP...)
                        theta = make_perm(decomp[1],decomp[2], decomp[3], index)
                        #println(theta)
                        #flush(stdout)
                        phi = psi^(-1)*theta^(-1)
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
            #end
		end
        #println("floop time = ")
        #println(string(time()-flooptime))
		return([Perm(Vector{Int}(sigma)), reduce(vcat,outlist)])
	end
end
