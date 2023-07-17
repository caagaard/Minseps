include("perm_utils.jl")
using Oscar

struct Hypermap
        sigma::Perm{Int}
        alpha::Perm{Int}
end

# Takes an integer i and permutations phi, psi.
# Returns [n,m] where n is the length of the orbit of i under phi and m is the length of the orbit of i
# under psi
function point_type(i::Int, phi::Perm{Int},psi::Perm{Int})
	return([[length(cycles(phi)[j] for j in 1:length(cycles(phi)) if i in cycles(phi)[j]][1],
		length(cycles(psi)[j]) for j in 1:length(cycles(psi)) if i in cycles(psi)[j]][1])
end

#Takes a pair of integers x and y and permutations phi and psi.  Returns true if there is a
# "color swapping isormophism" mapping x to y
function is_color_swap_iso(x::Int,y::Int,H::Hypermap)
	E = sum(conjclass(H.sigma))
	f_array = [0 for i in 1:E]
	# fill in f_array
	f_array[x] = y
	counter = 0
	while minimum(f_array)==0 && counter <= E
		for i in 1:E
			if f_array[i] != 0
				f_array[H.sigma[i]] = H.alpha[f_array[i]]
				f_array[H.alpha[i]] = H.sigma[f_array[i]]
			end
		end
		counter = counter +1
	end
	if sort(f_array) != [1:E;]
		return(0)
	end
	f = Perm(f_array)
		if f*H.sigma == H.alpha && f*H.alpha == H.sigma*f
			return(1)
		end
	return(0)
end

# Determins if a hypermap is "self-color-dual"
function is_self_color_dual(hypermap::Hypermap)
	E = sum(conjclass(H.sigma))
	# First easy check is is cycle structures match
	if conjclass(H.sigma) != conjclass(H.alpha)
		return(0)
	end
	# Now we check for isomorphism
	# Possible destinations for 1: set n = length of the cycle in sigma containing 1, m = lengt of 
	# cycle in alpha containing 1.  Destinations are all other elements with type (m,n)
	one_targets = [y for y in 2:E if point_type(y,H.sigma, H.alpha) == point_type(1, H.alpha, H.sigma)]
	for y in one_targets
		if(is_color_swap_iso(1,y,H)) ==1
			return(1)
		end
	end
	return(0)
end

# Does this method work? Verify and speed compare to other
#function is_color_swap_dual(x::Perm, y::Perm)
#	isdual =0
#	if conjclass(x) != conj(class(y)
#		return(0)
#	end
#	E = sum(conjclass(x))
#	t = solve_conjugation(x,y,E)
#	S = symmetric_group(E)
#	sigma = perm(S, [x[i] for i in 1:E])
#	alpha = perm(S, y[i] for i in 1:E])
#	phi = perm(S, [t[i] for i in 1:E])
#	H = centralizer(S,sigma)[1]
#	for h in H
#		if (phi^(-1)*h^(-1)*alpha*h*phi == alpha
#			isdual = 1
#			return(isdual)
#		end
#	end
#	return(isdual)
#end

function get_dual_map(psi::Perm{Int}, phi::Perm{Int}, n::Int)
        phitemp = [n+phi[i] for i in Int8(1):n]
        theta = Perm(vcat([psi[i] for i in Int8(1):n],phitemp))
        alpha = Perm(vcat([i+n for i in Int8(1):n],[i for i in Int8(1):n]))
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

function graph_from_embedding(nu::Perm{Int}, e::Int)
    v= length(cycles(nu))
    Adj = Matrix{Int64}(undef, v, v)
    for (i, vert) in enumerate(cycles(nu))
        matchers = flipVert(vert, e)
        for (j, altvert) in enumerate(cycles(nu))
            s = length(intersect(altvert, matchers))
            #print(s)
	    if i ==j 
		Adj[i,j] = div(s,2)
            else
                Adj[i,j] = s
	    end
        end
    end
    return(Adj)
end 

