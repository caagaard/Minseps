include("perm_utils.jl")
using Oscar

# Takes an integer i and permutations phi, psi.
# Returns [n,m] where n is the length of the orbit of i under phi and m is the length of the orbit of i
# under psi
function point_type(i::Int, phi::Perm{Int},psi::Perm{Int})
	return([[length(cycles(phi)[j]) for j in 1:length(cycles(phi)) if i in cycles(phi)[j]][1],
		[length(cycles(psi)[j]) for j in 1:length(cycles(psi)) if i in cycles(psi)[j]][1]])
end

#Takes a pair of integers x and y and permutations phi and psi.  Returns true if there is a
# "color swapping isormophism" mapping x to y
function is_color_swap_iso(x::Int,y::Int,sigma, alpha)
	E = sum(conjclass(sigma))
	f_array = [0 for i in 1:E]
	# fill in f_array
	f_array[x] = y
	counter = 0
	while minimum(f_array)==0 && counter <= E
		for i in 1:E
			if f_array[i] != 0
				f_array[sigma[i]] = alpha[f_array[i]]
				f_array[alpha[i]] = sigma[f_array[i]]
			end
		end
		counter = counter +1
	end
	if sort(f_array) != [1:E;]
		return(0)
	end
	f = Perm(f_array)
		if f*sigma == alpha*f && f*alpha == sigma*f
			return(1)
		end
	return(0)
end

# Determines if a hypermap is "self-color-dual"
function is_self_color_dual(sigma::Perm, alpha::Perm)
	E = sum(conjclass(sigma))
	# First easy check is is cycle structures match
	if conjclass(sigma) != conjclass(alpha)
		return(0)
	end
	# Now we check for isomorphism
	# Possible destinations for 1: set n = length of the cycle in sigma containing 1, m = lengt of 
	# cycle in alpha containing 1.  Destinations are all other elements with type (m,n)
	one_targets = [y for y in 1:E if point_type(y,sigma, alpha) == point_type(1, alpha, sigma)]
	for y in one_targets
		if(is_color_swap_iso(1,y,sigma, alpha)) ==1
			return(1)
		end
	end
	return(0)
end

# Takes a hypermap as input and returns a combinatorial map correspond to the dual hypermap
function get_dual_map(psi::Perm{Int}, phi::Perm{Int}, n::Int)
        phitemp = [n+phi[i] for i in Int(1):n]
        theta = Perm(vcat([psi[i] for i in Int(1):n],phitemp))
        alpha = Perm(vcat([i+n for i in Int(1):n],[i for i in Int(1):n]))
        sigma = theta^(-1)*alpha
        return([sigma,alpha])
end

function flipVert(vert, e::Int)
    flipped = []
    for i in vert
        if i >e
            j=i-e
        else
            j = i+e
        end
        push!(flipped, j)
    end
    return(flipped)
end 

# Takes a permutation nu in S_{2e} as input
# Constructs a graph homeomorphic to the underlying graph of the
# combinatorial map (nu, (1 e+1)(2 e+2)... (e 2e))
# Subdivides any duplicate edges to obtain a simple graph (with possible loops)
# This is because built in graph isomorphism method doesn't accept multigraphs
function graph_from_embedding(nu::Perm{Int}, e::Int)
    G = Graphs.Graph(length(cycles(nu)))
    for i in 1:length(cycles(nu))
	for x in cycles(nu)[i]
            if x < e+1
		for j in 1:length(cycles(nu))
		    if x+e in cycles(nu)[j]
		        if Graphs.has_edge(G, i, j)
			    if i == j
				Graphs.add_vertex!(G)
				Graphs.add_vertex!(G)
				Graphs.add_edge!(G, i, Graphs.nv(G)-1)
				Graphs.add_edge!(G, i, Graphs.nv(G))
				Graphs.add_edge!(G, Graphs.nv(G)-1, Graphs.nv(G))
				break
			    else
                            	Graphs.add_vertex!(G)
                            	Graphs.add_edge!(G, i, Graphs.nv(G))
                            	Graphs.add_edge!(G, j, Graphs.nv(G))
                            	break
			    end
                        else
                            Graphs.add_edge!(G,i,j)
                            break
                        end
                    end
                end
            end
        end
    end
    return(G)
end

