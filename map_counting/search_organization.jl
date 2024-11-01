include("hypermap_search.jl")
include("hypermap_utils.jl")
include("perm_utils.jl")
using Graphs
using Combinatorics
using DataStructures

# Returns the order of the conjugacy class in the symmetric group of a partition
#function conj_class_size(part::Vector{Int})
#	n = sum(part)
#	l = length(part)
	# The size of a conjugacy class is equal to the index of the centralizer
#	# Order of the centralizer is the product of the cycle lengths times the number of ways to mix them up.
#	numb_k_cycles = counter(part)
#	C_ord = prod(part)*prod([factorial(big(numb_k_cycles[k])) for k in 1:n])
#	return(factorial(big(n))/C_ord)
#end

#  First step: given g,ghat, E, find all valid triples lambda_1,lambda_2,lambda_3 such that they are minsep candidates.  Note at this point choice of one lambda doesn't affect the others, so just make 3 lists.
function get_class_candidates(g::Int, ghat::Int, E::Int, i::Int)
	# recall that the number of vertices is 2+g-ghat, and i is the number of black vertices.
	j = 2+g-ghat -i
	sigma_candidates = Combinatorics.partitions(E, i)
	alpha_candidates = Combinatorics.partitions(E, j)
	# F, the number of faces (cycles of phi) is given by 
	F = E - g - ghat
	# Candidate partitions for theta have no 1-cycles, so we construct partitions of E-F, and then add one to each block.
	phi_candidates = [x + ones(Int, F) for x in collect(Combinatorics.partitions((E-F),F))]
	candidates = [sigma_candidates, alpha_candidates, phi_candidates]
	# Now we need to compute some sizes
	n_sigma = sum([conj_class_size(part) for part in sigma_candidates])
	n_alpha = sum([conj_class_size(part) for part in alpha_candidates])
	n_phi = sum([conj_class_size(part) for part in phi_candidates])
	if n_sigma*length(alpha_candidates) >= n_alpha*length(sigma_candidates)
		if n_sigma >= n_phi
			theta_flag = 1
		else
			theta_flag = 0
		end
		#return(theta_flag, psi_candidates, theta_candidates, j, 0)
        return(theta_flag, sigma_candidates, phi_candidates, alpha_candidates)
	else
		if n_sigma >= n_phi
			theta_flag =1
		else
			theta_flag =0
		end
		#return(theta_flag, phi_candidates, theta_candidates, i, 1)
        # was the 5th value in this tuple ever used?
        return(theta_flag, alpha_candidates, phi_candidates, sigma_candidates)
	end
end

# Constructs a `default' permutation with a given cycle type
function make_default_perm(in_partition::Vector{Int})
	defperm = Vector{Vector{Int}}()
	Counter = 1
	for i in 1:length(in_partition)
		push!(defperm, [Counter:Counter+in_partition[i]-1;])
		Counter = Counter+in_partition[i]
	end
	return(defperm)
end

function get_ghat_minseps(g::Int, ghat::Int)
	big_ghat_minseps = Vector{Vector{Vector{Int}}}[]
	for E in (g+ghat+1):(2*(g+ghat))
		x = time()
		push!(big_ghat_minseps, get_ghat_minseps_edges(g, ghat, E))
		println(string(E))
		#flush(stdout)
        #println("E edges time =")
		#println(string(time()-x))
		#flush(stdout)
	end
    #println("about to reduce")
    #flush(stdout)
	ghat_minseps = reduce(vcat, big_ghat_minseps)
	return(ghat_minseps)
end

# Find minimal separating ribbon graphs with ribbon graph genus ghat and e edges
# Specialized version of get_ghat_minseps to be more easily split up for long calculations on multiple nodes
function get_ghat_minseps_edges(g::Int, ghat::Int, E::Int)
	ghat_minseps_E = Vector{Vector{Int}}[]
	needed_vertices = 2+g-ghat
	for i in 1:div(needed_vertices,2)
		conj_class_nums = get_class_candidates(g, ghat, E, i)
		sigma_choices = conj_class_nums[2] #Combinatorics.partitions(E, i)
		if conj_class_nums[1] ==1
			phi_choices = conj_class_nums[3]
			for sigma_choice in sigma_choices
				sigma = make_default_perm(sigma_choice)
                n_alpha_cycles = 2+g-ghat-length(sigma_choice)
				for phi_choice in phi_choices
                    if g-ghat >1
                        #need to fix conj_class_nums[4] here after change to conj_class_nums
						append!(ghat_minseps_E, [x for x in get_phi_candidates_v1(E,phi_choice, ghat, sigma, n_alpha_cycles,1)])
                    else
                        # need same fix here
                    	append!(ghat_minseps_E, [x for x in get_phi_candidates_v1(E,phi_choice, ghat, sigma, n_alpha_cycles,1,1)])
                    end
					#outs = reduce(vcat,tempouts)
					# Need to make sigma into an actual Perm{Int} object now:
					#S = symmetric_group(E)
					#append!(ghat_minseps_E, [x for x in outs])
				end
			end
		else
			#phi_cycles = conj_class_nums[4]
            alpha_choices = conj_class_nums[4]
            n_phi_cycles = E-g-ghat
			for sigma_choice in sigma_choices
			    sigma = make_default_perm(sigma_choice)
                for alpha_choice in alpha_choices
                    append!(ghat_minseps_E, [x for x in get_alphas_dist(E,ghat, sigma, alpha_choice, n_phi_cycles)])
				end
				#outs = get_phi_candidates_v2(E, phi_cycles, ghat, sigma)
				#append!(ghat_minseps_E, [x for x in outs])
            end
		end
	end
	#ghat_minseps_E = reduce(vcat, ghat_minseps_E)
    println("Length of minseps for ghat=$ghat and E=$E")
    println(string(length(ghat_minseps_E)))
    flush(stdout)
	if g-ghat > 1
		return([x for x in ghat_minseps_E if is_transitive_pair([Perm(x[1]), Perm(x[2])])])
	else
		return(ghat_minseps_E)
	end
	flush(stdout)
end


function generate_minseps_genus(g::Int)
	total_minseps = []
	for ghat in 0:g
		println("g, ghat = ")
		print(string(g))
		println(string(ghat))
		flush(stdout)
		glist = get_ghat_minseps(g,ghat)
		push!(total_minseps, glist)
	end
	return(reduce(vcat,total_minseps))
end

function count_embeds(hypermap_list::Vector{Vector{Vector{Int}}})
    tempcount = length(hypermap_list)
    self_color_counts = zeros(Threads.nthreads())
    Threads.@threads for hypermap in hypermap_list
        x = Perm(hypermap[1])
        y = Perm(hypermap[2])
	    if length(cycles(x)) == length(cycles(y))
		    if is_self_color_dual(x, y) == 0
			    #tempcount= tempcount -0.5
                self_color_counts[Threads.threadid()] += 1
		    end
	    end
    end
    tempcount = tempcount - 0.5*(sum(self_color_counts))
    return(tempcount)
end

function dual_list_to_minseps(dual_list::Vector{Vector{Vector{Int}}})
	minseps_list = [get_dual_map(Perm(dual[1]),Perm(dual[2]), Int(sum(length(k) for k in cycles(Perm(dual[1]))))) for dual in dual_list]
	return(minseps_list)
end

function getIsoClasses(graphlist)
    isoclasses = Vector{SimpleGraph}() 
    for graph in graphlist
        isdupe = 0
        Threads.@threads for G in isoclasses
            if Graphs.Experimental.has_isomorph(graph, G)==1
                isdupe = 1
                break
            end
            #extra check to avoid instability with threads
            if isdupe ==1
                break
            end
        end
        if isdupe == 0
            push!(isoclasses, graph)
        end
    end
    return(isoclasses)
end

function minseps_list_to_graphs(minsep_list::Vector{Vector{Perm{Int}}}, g::Int)
        temp_graphs_list = [graph_from_embedding(ribbon_graph[1], length(cycles(ribbon_graph[2]))) for ribbon_graph in minsep_list]
        sorted_graphs_list =[Vector{SimpleGraph}() for e in 1:4*g]
        for ribbon_graph in minsep_list
		        E = length(cycles(ribbon_graph[2]))
	            push!(sorted_graphs_list[E], graph_from_embedding(ribbon_graph[1], E))
        end
        # Need to deal with edge countness
        final_graphs = [getIsoClasses(sorted_graphs_list[e]) for e in 1:length(sorted_graphs_list)]
        fgs = reduce(vcat, final_graphs)
        return(fgs)
end

#function hypermap_to_dual_graph(sigma::Vector{Vector{Int}}, alpha::Vector{Vector{Int}})
#    dual_map = get_dual_map(Perm(sigma), Perm(alpha), Int(sum(length(k) for k in cycles(Perm(sigma)))))
#    return(graph_from_embedding(dual_map[1], length(cycles(dual_map[2]))))
#end
