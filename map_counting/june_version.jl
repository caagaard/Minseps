using Oscar
using Combinatorics
using FLoops

function perm_components(part::Array{Array{Int8,1},1}, n::Int8)
        starters = [chunk[1] for chunk in part]
        perm_parts = [chunk[2:length(chunk)] for chunk in part]
        return(starters, perm_parts,n)
end
function perm_counter(s::Array{Int64,1})
        return([1:factorial(i-1) for i in s])
end

function make_perm(starters::Array{Int8,1}, perm_parts::Array{Array{Int8,1},1}, n::Int8, index::NTuple)
        phi = [Int8(1):n;]
        for i in 1:length(starters)
                tempperm = vcat(starters[i], nthperm(perm_parts[i], index[i]))
                for j in 1:(length(tempperm)-1)
                        phi[tempperm[j]] = tempperm[j+1]
                        phi[tempperm[end]] = tempperm[1]
                end
        end
        return(Perm(phi))
end

function get_phi_candidates_v2(n::Int64,k::Int64,g::Int64, psitemp::Vector{Vector{Int64}})
	S = symmetric_group(n)
	sigma = cperm(S, psitemp...)
	psi = Perm(Vector{Int8}(sigma))
	needed_verts = n-k -length(cycles(psi))+2 - 2*g
	#println(psi)
	H = centralizer(S,sigma)
	HH = [Perm(Vector{Int8}(x)) for x in H[1]]
        parts = [part for part in Combinatorics.partitions([Int8(1):Int8(n);],k)]
	outlist = []
	for i in 1:Threads.nthreads()
		push!(outlist, [])
	end
        for part in parts
                decomp = perm_components(part,Int8(n))
                PP = perm_counter([length(part[i]) for i in 1:length(part)])
                flooptime = time()
                @floop for index in Iterators.product(PP...)
                        phi = make_perm(decomp[1],decomp[2], decomp[3], index)
			            is_min =1
			            for g in HH
				                if vcat(cycles(g^(-1)*phi*g)...) < vcat(cycles(phi)...)
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
        return([Perm(Vector{Int8}(sigma)), reduce(vcat,outlist)])
end 

function get_phi_candidates_v1(n::Int64, part::Vector{Int64}, g::Int64, psitemp::Vector{Vector{Int64}}, n_phi_cycles::Int64)
	if sum(part) != n
		println("Invalid Partition")
		return([])
	else
		S = symmetric_group(n)
		sigma = cperm(S,psitemp...)
		psi = Perm(Vector{Int8}(sigma))
		H = centralizer(S,sigma)
		HH = [Perm(Vector{Int8}(x)) for x in H[1]]
		outlist = []
		for i in 1:Threads.nthreads()
			push!(outlist, [])
		end
		thetap = cperm(S, make_default_perm(part)...)
		CC = [Perm(Vector{Int8}(x)) for x in conjugacy_class(S, thetap)]
		cc_dict = Dict([(c,i) for (i,c) in enumerate(CC)])
        flooptime = time()
		@floop for theta in CC
			is_min = 1
			for g in HH
				if cc_dict[theta^g] < cc_dict[theta]
					is_min=0
					break
				end
			end
			if is_min ==1
				phi = psi^(-1)*theta^(-1)
				# Testing with ==1 for cellular case
				#if length(cycles(phi)) == 1
				if length(cycles(phi)) == n_phi_cycles
					push!(outlist[Threads.threadid()], phi)
				end
			end
		end
        #println("floop time = ")
        #println(string(time()-flooptime))
		return([Perm(Vector{Int8}(sigma)), reduce(vcat,outlist)])
	end
end
