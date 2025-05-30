include("hypermap_search.jl")
include("hypermap_utils.jl")
include("perm_utils.jl")
using Graphs
using Combinatorics
using DataStructures

#  First step: given g,ghat, E, find all valid triples lambda_1,lambda_2,lambda_3 such that 
#  they are minsep candidates.  Note at this point choice of one lambda doesn't affect the 
#  others, so just make 3 lists.
#  Additional input i is the number of black vertices to break down the search into 
#  more manageable sized parts
#  Returns a list of conjugacy class types to consider.
#  phi_flag indicates whether the 3rd or 4th output is the smaller set to search
#  second element returned is the set from which we will fix an element
#  This corresponds to being able to freely assign labels to one element of the hypermap
function get_class_candidates(g::Int, ghat::Int, E::Int, i::Int)
	# recall that the number of vertices is 2+g-ghat.  i is the number of black vertices.
	# j will be the number of white vertices
	# In internal use we always have j>= i.  If you use this method without satisfying
	# that criterion, be very careful in how you use this method
	j = 2+g-ghat -i
	# sigma_candidates is the number of cycle types sigma could have, etc
	sigma_candidates = Combinatorics.partitions(E, i)
	alpha_candidates = Combinatorics.partitions(E, j)
	# F, the number of faces (cycles of phi) is given by 
	F = E - g - ghat
	# Candidate partitions for theta have no 1-cycles, so we construct partitions of E-F, and then add one to each block.
	phi_candidates = [x + ones(Int, F) for x in collect(Combinatorics.partitions((E-F),F))]
	candidates = [sigma_candidates, alpha_candidates, phi_candidates]
	# Now we need to compute some sizes n_sigma is the number of permutations which could be sigma, etc
	n_sigma = sum([conj_class_size(part) for part in sigma_candidates])
	n_alpha = sum([conj_class_size(part) for part in alpha_candidates])
	n_phi = sum([conj_class_size(part) for part in phi_candidates])
	# Check if it's faster to fix an element for each type for sigma and search all alpha candidates or vice-versa
	if n_sigma*length(alpha_candidates) >= n_alpha*length(sigma_candidates)
		# If faster to search phis than sigmas, set phi_flag=1
		if n_sigma >= n_phi
			phi_flag = 1
		else
			phi_flag = 0
		end
        	return(phi_flag, sigma_candidates, phi_candidates, alpha_candidates)
	else
		# If faster to search phis than alphas, set phi_flag=1
		if n_alpha >= n_phi
			phi_flag =1
		else
			phi_flag =0
		end
        	return(phi_flag, alpha_candidates, phi_candidates, sigma_candidates)
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

# Computes the a list of isomorphism classes of hypermaps which
# have genus ghat, E edges, and correspond to a ribbon graph in E_g
function get_ghat_minseps(g::Int, ghat::Int)
	big_ghat_minseps = Vector{Vector{Vector{Int}}}[]
	for E in (g+ghat+1):(2*(g+ghat))
		x = time()
		push!(big_ghat_minseps, get_ghat_minseps_edges(g, ghat, E))
		println(string(E))
		flush(stdout)
	end
	ghat_minseps = reduce(vcat, big_ghat_minseps)
	return(ghat_minseps)
end

# Find minimal separating ribbon graphs with ribbon graph genus ghat and e edges
# Specialized version of get_ghat_minseps to be more easily split up for long calculations on multiple nodes
function get_ghat_minseps_edges(g::Int, ghat::Int, E::Int)
	ghat_minseps_E = Vector{Vector{Int}}[]
	# The dual of a ribbon graph in E_g with ribbon graph genus ghat has 2+g-ghat vertices
	needed_vertices = 2+g-ghat
	# To avoid double counting we require the number of black vertices is at most the number of white vertices
	for i in 1:div(needed_vertices,2)
		conj_class_nums = get_class_candidates(g, ghat, E, i)
		sigma_choices = conj_class_nums[2] 
		# If conj_class_nums[1] =1, we search possible phis
		if conj_class_nums[1] ==1
			phi_choices = conj_class_nums[3]
			for sigma_choice in sigma_choices
				sigma = make_default_perm(sigma_choice)
				# Set number of cycles alpha must have to form a hypermap of correct genus
 				n_alpha_cycles = 2+g-ghat-length(sigma_choice)
				for phi_choice in phi_choices
					# Memory usage is not an issue when g-ghat>1, so we use a simpler search computing characters to preallocate space
 					if g-ghat >1
                        			#need to fix conj_class_nums[4] here after change to conj_class_nums
						append!(ghat_minseps_E, [x for x in get_phi_candidates_v1(E,phi_choice, ghat, sigma, n_alpha_cycles,1)])
					# When g-ghat <=1 some cases have very high memory demands so we use a search method that computes group characters
					# and uses a formula of Frobenius to preallocate space
 					else
 						append!(ghat_minseps_E, [x for x in get_phi_candidates_v1(E,phi_choice, ghat, sigma, n_alpha_cycles,1,1)])
 					end
				end
			end
		else
		# If conj_class_nums[1] != 1, we search possible alphas
			alpha_choices = conj_class_nums[4]
			# Set number of cycles phi must have to form a hypermap of correct genus
 			n_phi_cycles = E-g-ghat
			for sigma_choice in sigma_choices
				sigma = make_default_perm(sigma_choice)
 				for alpha_choice in alpha_choices
 					# Want to avoid finding both colorings of the same object
					# This is partially handled by setting i <= j. To avoid duplicate search when i=j
					# we don't search if the cycle type of alpha comes before that of sigma in lex order
 					if  2*i != 2+g-ghat || sigma_choice <= alpha_choice
 						append!(ghat_minseps_E, [x for x in get_alphas_dist(E,ghat, sigma, alpha_choice, n_phi_cycles)])
 					end
				end
 			end
		end
	end
 	println("Length of minseps for ghat=$ghat and E=$E")
 	println(string(length(ghat_minseps_E)))
	flush(stdout)

	# If g-ghat >1 we need to verify that <sigma, alpha> acts transitvely on {1,...,n} for 
	# the pair to form a hypermap
	if g-ghat > 1
		return([x for x in ghat_minseps_E if is_transitive_pair([Perm(x[1]), Perm(x[2])])])

	# If g-ghat <= 1 then by our construction sigma is an n cycle and the action is automatically transitive
	else
		return(ghat_minseps_E)
	end
end

# Generates a list of isomorphism classes of hypermaps
# which are dual to a ribbon graph in E_g
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

# Generates a list of isomorphism classes of hypermaps
# which have E edges and are dual to a ribbon graph in E_g
function generate_minseps_genus(g::Int, E::Int)
    total_minseps = []
    for ghat in maximum([0,(div(E,2)-g)]): minimum([g,(E-g-1)])
        println("g, ghat = ")
        println(string(g))
        println(string(ghat))
        flush(stdout)
        push!(total_minseps, get_ghat_minseps_edges(g, ghat, E))
    end
    return(reduce(vcat, total_minseps))
end

# Input is a list of hypermaps
# Output is a list of combinatorial maps, each one the combinatorial
# map corresponding to the dual hypermap of one of the inputs
function dual_list_to_minseps(dual_list::Vector{Vector{Vector{Int}}})
	minseps_list = [get_dual_map(Perm(dual[1]),Perm(dual[2]), Int(sum(length(k) for k in cycles(Perm(dual[1]))))) for dual in dual_list]
	return(minseps_list)
end

# Takes a list of graphs as input
# Returns a list with a single representative of each isomorphism 
# class present in input list
# Note: has_isomorph requires graphs are processed to avoid multiedges
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

# Takes a list of combinatorial maps as input
# Constructs a graph homeomorphic to the underlying graph of each map
# Returns a vector, whose i-th entry contains a single representative
# of each isomorphism class of graph underlying the input maps with 
# exactly i edges
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
