include("hypermap_search.jl")
include("hypermap_utils.jl")
include("perm_utils.jl")
using Graphs
using Combinatorics
using DataStructures

# Returns the order of the conjugacy class in the symmetric group of a partition
function conj_class_size(part::Vector{Int})
	n = sum(part)
	l = length(part)
	# The size of a conjugacy class is equal to the index of the centralizer
	# Order of the centralizer is the product of the cycle lengths times the number of ways to mix them up.
	numb_k_cycles = counter(part)
	C_ord = prod(part)*prod([factorial(big(numb_k_cycles[k])) for k in 1:n])
	return(factorial(big(n))/C_ord)
end

#  First step: given g,ghat, E, find all valid triples lambda_1,lambda_2,lambda_3 such that they are minsep candidates.  Note at this point choice of one lambda doesn't affect the others, so just make 3 lists.
function get_class_candidates(g::Int, ghat::Int, E::Int, i::Int)
	# recall that the number of vertices is 2+g-ghat, and i is the number of black vertices.
	j = 2+g-ghat -i
	psi_candidates = Combinatorics.partitions(E, i)
	phi_candidates = Combinatorics.partitions(E, j)
	# Now we determine the number of faces, F
	F = E - g - ghat
	# Candidate partitions for theta have no 1-cycles, so we construct partitions of E-F, and then add one to each block.
	theta_candidates = [x + ones(Int, F) for x in collect(Combinatorics.partitions((E-F),F))]
	candidates = [psi_candidates, phi_candidates, theta_candidates]
	# Now we need to compute some sizes
	n_psi = sum([conj_class_size(part) for part in psi_candidates])
	n_phi = sum([conj_class_size(part) for part in phi_candidates])
	n_theta = sum([conj_class_size(part) for part in theta_candidates])
	if n_psi*length(phi_candidates) >= n_phi*length(psi_candidates)
		if n_phi > n_theta
			theta_flag = 1
		else
			theta_flag = 0
		end
		return(theta_flag, psi_candidates, theta_candidates, j, 0)
	else
		if n_psi > n_theta
			theta_flag =1
		else
			theta_flag =0
		end
		return(theta_flag, phi_candidates, theta_candidates, i, 1)
	end
end

function make_default_perm(in_partition::Vector{Int})
	defperm = Vector{Vector{Int}}()
	Counter = 1
	for i in 1:length(in_partition)
		push!(defperm, [Counter:Counter+in_partition[i]-1;])
		Counter = Counter+in_partition[i]
	end
	#println(typeof(defperm))
	flush(stdout)
	return(defperm)
end

function get_ghat_minseps(g::Int, ghat::Int)
	big_ghat_minseps = Vector{Vector{Perm{Int}}}[]
	for E in (g+ghat+1):(2*(g+ghat))
		x = time()
		push!(big_ghat_minseps, get_ghat_minseps_edges(g, ghat, E))
		println(string(E))
		#flush(stdout)
        	println("E edges time =")
		println(string(time()-x))
		flush(stdout)
	end
	ghat_minseps = reduce(vcat, big_ghat_minseps)
	return(ghat_minseps)
end

# Find minimal separating ribbon graphs with ribbon graph genus ghat and e edges
# Specialized version of get_ghat_minseps to be more easily split up for long calculations on multiple nodes
function get_ghat_minseps_edges(g::Int, ghat::Int, E::Int)
	ghat_minseps_E = Vector{Perm{Int}}[]
	needed_vertices = 2+g-ghat
	for i in 1:div(needed_vertices,2)
		conj_class_nums = get_class_candidates(g, ghat, E, i)
		psi_choices = conj_class_nums[2] #Combinatorics.partitions(E, i)
		if conj_class_nums[1] ==1
			theta_choices = conj_class_nums[3]
			for psi_choice in psi_choices
				psi = make_default_perm(psi_choice)
				for theta_choice in theta_choices
					append!(ghat_minseps_E, [x for x in get_phi_candidates_thread(E,theta_choice, ghat, psi, (conj_class_nums[4]))])
				end
				#outs = reduce(vcat,tempouts)
				# Need to make psi into an actual Perm{Int} object now:
<<<<<<< HEAD
				S = symmetric_group(E)
				push!(ghat_minseps_E, [[Perm(Vector{Int}(cperm(S,psi...))), phi] for phi in outs])
=======
				#S = symmetric_group(E)
				#append!(ghat_minseps_E, [x for x in outs])
>>>>>>> df4858a58ea61eb55587a6052252a905d086774a
			end
		else
				phi_cycles = conj_class_nums[4]
				for psi_choice in psi_choices
					psi = make_default_perm(psi_choice)
					outs = get_phi_candidates_v2(E, phi_cycles, ghat, psi)
<<<<<<< HEAD
					push!(ghat_minseps_E, [[outs[1], phi] for phi in outs[2]])
				end
		end
	end
    ghat_minseps_E = reduce(vcat, ghat_minseps_E)
    flush(stdout)
    if g-ghat > 1
        connected_minseps = [x for x in ghat_minseps_E if is_transitive_pair(x)]
        return(connected_minseps)
    else
        return(ghat_minseps_E)
    end
	#Counter = 0
	#for psi_choice in ghat_minseps_E
	#	if length(psi_choice) == 2
	#		Counter = Counter+length(psi_choice[2])
	#	else
	#		println("unexpected array length")
	#		println(length(psi_choice))
	#	end
	#end
	#flush(stdout)
	#return(ghat_minseps_E)
    #return(connected_minseps)
=======
					append!(ghat_minseps_E, [x for x in outs])
				end
		end
	end
	#ghat_minseps_E = reduce(vcat, ghat_minseps_E)
	if g-ghat > 1
		return([x for x in ghat_minseps_E if is_transitive_pair(x)])
	else
		return(ghat_minseps_E)
	end
	flush(stdout)
>>>>>>> df4858a58ea61eb55587a6052252a905d086774a
end


function generate_minseps_genus(g::Int)
	total_minseps = []
	for ghat in 0:g
		println("g, ghat = ")
		print(string(g))
		println(string(ghat))
		flush(stdout)
		glist = get_ghat_minseps(g,ghat)
<<<<<<< HEAD
		#biggerlist = []
		#for x in glist
			#templist = [[x[1], i] for i in x[2]]
            #println(string(length(templist)))
            #flush(stdout)
			#push!(biggerlist, x)
		#end
		#ghypermaps = reduce(vcat,glist)
		#ghypermaps = [permpair for permpair in gbiglist if is_transitive_pair(permpair)]
=======
>>>>>>> df4858a58ea61eb55587a6052252a905d086774a
		push!(total_minseps, glist)
	end
	return(reduce(vcat,total_minseps))
end

function count_embeds(hypermap_list::Vector{Vector{Perm{Int}}})
    tempcount = length(hypermap_list)
    for hypermap in hypermap_list
	if length(cycles(hypermap[1])) == length(cycles(hypermap[2]))
		if is_self_color_dual(hypermap[1], hypermap[2]) == 0
			tempcount= tempcount -0.5
		end
	end
    end
    return(tempcount)
end
function dual_list_to_minseps(dual_list::Vector{Vector{Perm{Int}}})
	minseps_list = [get_dual_map(dual[1],dual[2], Int(sum(length(k) for k in cycles(dual[1])))) for dual in dual_list]
	return(minseps_list)
end

function getIsoClasses(graphlist)
        isoclasses = Vector{SimpleGraph}() 
        for graph in graphlist
                isdupe = 0
                for G in isoclasses
                        if Graphs.Experimental.has_isomorph(graph, G)==1
                                isdupe = 1
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
