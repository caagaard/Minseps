using Combinatorics
using DataStructures
using Oscar

# Returns the conjugacy class of a permutation as a list of cycle lengths
function conjclass(x::Perm{Int})
	return(sort([length(i) for i in cycles(x)]))
end

function conj_class_size(part::Vector{Int})
	n = sum(part)
	l = length(part)
	# The size of a the conjugacy class in S_n is equal to the index of the centralizer subgroup
	# of any element of the conjugacy class
	# The order of the centralizer is the product of the cycle lengths times the number of ways to 
	# mix up the cycles of the same length
	numb_k_cycles = counter(part)
	C_ord = prod(part)*prod([factorial(big(numb_k_cycles[k])) for k in 1:n])
	return(factorial(big(n))/C_ord)
end

# If two permutations x and y are conjugate, returns a permutation t such that x^t = y.  Otherwise 
# returns 0.
function solve_conjugation(x::Perm{Int}, y::Perm{Int}, n::Int)
	x_cyc_type = sort([length(cycles(x)[i]) for i in 1:length(cycles(x))])
	y_cyc_type = sort([length(cycles(y)[i]) for i in 1:length(cycles(y))])
	if x_cyc_type != y_cyc_type
		return(0)
	else
		VV = vcat(sort([cycles(x)[i] for i in 1:length(cycles(x))], by=length)...)
		WW = vcat(sort([cycles(y)[i] for i in 1:length(cycles(y))], by=length)...)
		tau = zeros(Int64, n)
		for i in 1:n
			tau[VV[i]] = WW[i]
		end
		t = Perm(tau)
		if t^-1 *x*t != y
			println("OH NO!")
		end
		return(t)
	end
end

# Makes a "default" permutation with given cycle type.  Returns it as a list of integers
function make_default_perm(in_partition::Vector{Int})
	defperm = Vector{Vector{Int64}}()
	counter = 1
	for i in 1:length(in_partition)
		push!(defperm, [counter:counter+in_partition[i]-1;])
		counter = counter+in_partition[i]
	end
	return(defperm)
end

# Takes a pair of permutations and returns whether they generate a transitive permutation group
# Not sure about typing on the input
function is_transitive_pair(permpair::Vector{Perm{Int}})
	n = sum(conjclass(permpair[1]))
	S = SymmetricGroup(n)
	sigma = cperm(S, [cycles(permpair[1])[i] for i in 1:length(cycles(permpair[1]))])
	alpha = cperm(S, [cycles(permpair[2])[i] for i in 1:length(cycles(permpair[2]))])
	# type adjust sigma and alph to be elements of S
	H = Oscar.permutation_group(n, [sigma, alpha])
	return(is_transitive(H))
end	

