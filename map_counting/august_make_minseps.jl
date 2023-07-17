include("june_version.jl")
include("hypermap_duals.jl")
#include("graphmethods.jl")
include("check_transitive.jl")
using Graphs
using Combinatorics
using DataStructures

# Returns the order of the conjugacy class in the symmetric group of a partition
function conj_class_size(part::Vector{Int64})
	n = sum(part)
	l = length(part)
	# The size of a conjugacy class is equal to the index of the centralizer
	# Order of the centralizer is the product of the cycle lengths times the number of ways to mix them up.
	numb_k_cycles = counter(part)
	C_ord = prod(part)*prod([factorial(big(numb_k_cycles[k])) for k in 1:n])
	return(factorial(big(n))/C_ord)
end

#  First step: given g,ghat, E, find all valid triples lambda_1,lambda_2,lambda_3 such that they are minsep candidates.  Note at this point choice of one lambda doesn't affect the others, so just make 3 lists.
function get_class_candidates(g::Int64, ghat::Int64, E::Int64, i::Int64)
	# recall that the number of vertices is 2+g-ghat, and i is the number of black vertices.
	j = 2+g-ghat -i
	psi_candidates = Combinatorics.partitions(E, i)
	phi_candidates = Combinatorics.partitions(E, j)
	# Now we determine the number of faces, F
	F = E - g - ghat
	# Candidate partitions for theta have no 1-cycles, so we construct partitions of E-F, and then add one to each block.
	theta_candidates = [x + ones(Int64, F) for x in collect(Combinatorics.partitions((E-F),F))]
	candidates = [psi_candidates, phi_candidates, theta_candidates]
	# Now we need to compute some sizes
	n_psi = sum([conj_class_size(part) for part in psi_candidates])
	n_phi = sum([conj_class_size(part) for part in phi_candidates])
	n_theta = sum([conj_class_size(part) for part in theta_candidates])
	if n_phi > n_theta
		theta_flag = 1
	else
		theta_flag = 0
	end
	return(theta_flag, psi_candidates, theta_candidates, j)
end

function make_default_perm(in_partition::Vector{Int64})
	defperm = Vector{Vector{Int64}}()
	Counter = 1
	for i in 1:length(in_partition)
		push!(defperm, [Counter:Counter+in_partition[i]-1;])
		Counter = Counter+in_partition[i]
	end
	return(defperm)
end

function get_ghat_minseps(g::Int64, ghat::Int64)
	ghat_minseps = []
	for E in (g+ghat+1):(2*(g+ghat))
	#for E in (2*(g+ghat)):(2*(g+ghat))
		x = time()
		push!(ghat_minseps, get_ghat_minseps_edges(g, ghat, E))
		println(string(E))
		flush(stdout)
        println("E edges time =")
		println(string(time()-x))
		flush(stdout)
	end
	ghat_minseps = reduce(vcat, ghat_minseps)
	return(ghat_minseps)
end

# Find minimal separating ribbon graphs with ribbon graph genus ghat and e edges
# Specialized version of get_ghat_minseps to be more easily split up for long calculations on multiple nodes
function get_ghat_minseps_edges(g::Int64, ghat::Int64, E::Int64)
	ghat_minseps_E = []
	needed_vertices = 2+g-ghat
	for i in 1:div(needed_vertices,2)
		conj_class_nums = get_class_candidates(g, ghat, E, i)
		# Remove later, this line for test now:
		#conj_class_nums[1]=0
		psi_choices = Combinatorics.partitions(E, i)
		if conj_class_nums[1] ==1
			theta_choices = conj_class_nums[3]
			for psi_choice in psi_choices
				psi = make_default_perm(psi_choice)
				tempouts = []
				for theta_choice in theta_choices
					push!(tempouts, get_phi_candidates_v1(E,theta_choice, ghat, psi, (needed_vertices-i))[2])
				end
				outs = reduce(vcat,tempouts)
				# Need to make psi into an actual Perm{Int8} object now:
				S = symmetric_group(E)
				push!(ghat_minseps_E, [Perm(Vector{Int8}(cperm(S,psi...))), outs])
			end
		else
			phi_cycles = conj_class_nums[4]
			for psi_choice in psi_choices
				psi = make_default_perm(psi_choice)
				outs = get_phi_candidates_v2(E, phi_cycles, ghat, psi)
				push!(ghat_minseps_E, outs)
			end
		end
	end
	Counter = 0
	for psi_choice in ghat_minseps_E
		if length(psi_choice) == 2
			Counter = Counter+length(psi_choice[2])
		else
			println("unexpected array length")
			println(length(psi_choice))
		end
	end
	#println(string(counter))
	#println(string(ghat_minseps_E))
	flush(stdout)
	return(ghat_minseps_E)
#	biggerlist = []
#	for x in ghat_minseps_E
#		templist = [[x[1], i] for i in x[2]]
#		push!(biggerlist, templist)
#	end
#	gbiglist = vcat(biggerlist...)
#	ghypermaps_E = [permpair for permpair in gbiglist if is_trans(permpair)]
#	return(ghypermaps_E)
end


function generate_minseps_genus(g::Int64)
	total_minseps = []
	for ghat in 0:g
#	for ghat in [g]
		println("g, ghat = ")
		print(string(g))
		println(string(ghat))
		flush(stdout)
		glist = get_ghat_minseps(g,ghat)
		biggerlist = []
		for x in glist
			templist = [[x[1], i] for i in x[2]]
			push!(biggerlist, templist)
		end
		gbiglist = reduce(vcat,biggerlist)
		ghypermaps = [permpair for permpair in gbiglist if is_trans(permpair)]
		push!(total_minseps, ghypermaps)
	end
	return(reduce(vcat,total_minseps))
end

# Returns the conjugacy class of a permutation as a list of cycle lengths
function conjclass(x::Perm{Int8})
	return(sort([length(i) for i in cycles(x)]))
end

function solve_conjugation(x::Perm{Int8}, y::Perm{Int8}, n)
           VV = reduce(vcat,sort([cycles(x)[i] for i in 1:length(cycles(x))], by=length))
           WW = reduce(vcat,sort([cycles(y)[i] for i in 1:length(cycles(x))], by=length))
           x_cyc_type = sort([length(cycles(x)[i]) for i in 1:length(cycles(x))])
           y_cyc_type = sort([length(cycles(y)[i]) for i in 1:length(cycles(y))])
           if x_cyc_type != y_cyc_type
               return(0)
           else
               tau = zeros(Int64, n)
               for i in 1:n
                   tau[VV[i]] =WW[i]
               end
               t = Perm(tau)
               if t^-1 * x * t != y
                   println("OH NO!")
               end
               return(t)
           end
end

# Takes a hypermap as input and determines if color-swapping is an isomorphism
function is_self_dual(x::Perm{Int8}, y::Perm{Int8})
	isdual = 0
	if conjclass(x) != conjclass(y)
		return(0)
	end
	E = sum(conjclass(x))
	t = solve_conjugation(x,y,E)
	S = symmetric_group(E)
	psi = perm(S, [x[i] for i in 1:E])
	phi = perm(S, [y[i] for i in 1:E])
	rho = perm(S, [t[i] for i in 1:E])
	H = centralizer(S, psi)[1]
	for h in H
		if (rho^(-1)*h^(-1)*phi*h*rho) == psi
			isdual = 1
			return(isdual)
		end
	end
	return(isdual)
end


#function self_duals_check(minseps_list)
#	dif_conj = 0
#	over_counter = 0
#	self_dual_counter = 0
#	for minsep in minseps_list
#		psi = minsep[1]
#		phi = minsep[2]
#		if conjclass(psi) == conjclass(phi)
#			if is_self_dual(psi, phi)==1
#				self_dual_counter = self_dual_counter +1
#			else
#				over_counter = over_counter +1
#			end
#		else
#			dif_conj = dif_conj+1
#		end
#	end
#	println("self duals is ")
#	println(string(self_dual_counter))
#	println("dif conj classes is ")
#	println(string(dif_conj))
#	flush(stdout)
#	return(over_counter)
#end

function dual_list_to_minseps(dual_list::Vector{Vector{Perm{Int8}}})
	minseps_list = [get_dual_map(dual[1],dual[2], Int8(sum(length(k) for k in cycles(dual[1])))) for dual in dual_list]
#	for dual in dual_list
#		temp_minseps_list = [get_dual_map(dual[1], dual[2], Int8(sum(length(k) for k in cycles(dual[1])))) for i in 1:length(dual[2])]
#		minseps_list = vcat(minseps_list, temp_minseps_list)
#	end
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

function minseps_list_to_graphs(minsep_list::Vector{Vector{Perm{Int8}}}, g::Int)
        temp_graphs_list = [graph_from_embedding(ribbon_graph[1], length(cycles(ribbon_graph[2])),1) for ribbon_graph in minsep_list]
        sorted_graphs_list =[Vector{SimpleGraph}() for e in 1:4*g]
        println("Hello?")
        flush(stdout)
        for ribbon_graph in minsep_list
		        E = length(cycles(ribbon_graph[2]))
	            push!(sorted_graphs_list[E], graph_from_embedding(ribbon_graph[1], E,1))
        end
        # Need to deal with edge countness
        final_graphs = [getIsoClasses(sorted_graphs_list[e]) for e in 1:length(sorted_graphs_list)]
        fgs = reduce(vcat, final_graphs)
        return(fgs)
end
