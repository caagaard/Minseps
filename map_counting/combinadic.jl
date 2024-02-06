# Implements ranking unranking algorithms for k-combinations from sets of size n
# Assumes that
using SplittablesBase
using DataStructures

function combination_next(n::Int,k::Int, v::Vector{Int})
        j=k
        while v[j] == n-k+j
                j=j-1
                if j==0
                        return(v)
                end
        end
        v[j] = v[j]+1
        for i in (j+1):k
                v[i] = v[i-1]+1
        end
        return(v)
end

# length_vec should be the following vector: [sum([length(x) for x in v[i:end]]) for i in 1:length(v)]
# Precomputing this reduces allocations by about 50%
function combination_tuple_next(v::Vector{Vector{Int}}, length_vec::Vector{Int})
        for i in length(v):-1:1
                k = length(v[i])
                if v[i] != [1; length_vec[i]+2-length(v[i]):length_vec[i];]
                        #println(i)
                        v[i] = combination_next(length_vec[i],k, v[i])
                        #if i > 1 && i <length(v)
                         #       if (length(v[i]) == length(v[i-1]) && v[i][1] < v[i-1][1]) || (length(v[i]) ==length(v[i+1]) && v[i][1] > v[i+1][1])
                         #               return(combination_tuple_next(v, length_vec))
                         #       end
                        #end
                        return(v)
                end
                v[i] = [1:k;]
        end
        return(v)
end

function regular_combination_next!(v::Vector{Vector{Int}}, k::Int)
        n = length(v)*k
        length_vec = [n:-k:k;]
        for i in length(v):-1:1
                if v[i] != [1; length_vec[i]+2-k:length_vec[i];]
                        #println(i)
                        v[i] = combination_next(length_vec[i],k, v[i])
                        #if i > 1 && i <length(v)
                         #       if (length(v[i]) == length(v[i-1]) && v[i][1] < v[i-1][1]) || (length(v[i]) ==length(v[i+1]) && v[i][1] > v[i+1][1])
                         #               return(combination_tuple_next(v, length_vec))
                         #       end
                        #end
                        return(v)
                end
                v[i] = [1:k;]
        end
        return(v)
end

struct RegularCombinationIterator{T}
        data::T
        start::T
        l::Int
        r::Int
        Length::Int
end

function Base.iterate(rc::RegularCombinationIterator)
        rc.data, 1
end

function Base.iterate(rc::RegularCombinationIterator, state)
        regular_combination_next!(rc.data, rc.l)
        state += 1
        if state > length(rc)
                return nothing
        end
        rc.data, state
end

Base.length(rc::RegularCombinationIterator) = rc.Length
Base.eltype(rc::RegularCombinationIterator) = Vector{Vector{Int}}

function makehalve(RC::RegularCombinationIterator)#{Vector{Vector{Int64}}})
        LL = div(RC.Length,2)
        left = RegularCombinationIterator(RC.data, RC.start, RC.l, RC.r, LL)
        right = RegularCombinationIterator(unrank_reg_combo(rank_reg_combo(RC.start, RC.l, RC.r)+LL, RC.l*RC.r, RC.l),
                unrank_reg_combo(rank_reg_combo(RC.start, RC.l, RC.r)+LL, RC.l*RC.r, RC.l), RC.l, RC.r, RC.Length - LL)
        return(left, right)
end

SplittablesBase.halve(RC::RegularCombinationIterator) = makehalve(RC)

function makeRCI(l::Int, r::Int)
        L = div(factorial(r*l), (factorial(r)*factorial(l)^r))
        return(RegularCombinationIterator([[1:l;] for i in 1:r], [[1:l;] for i in 1:r], l, r, L))
end

function unrank_reg_combo(r::Int, n::Int, k::Int)
        c = div(n,k)
        if c*k != n
                println("error k does not divide n")
                return([0])
        elseif n <=k
                #println("n=0")
                return([[1:k;]])
        else
                Q = div(factorial(n-k), factorial(c-1)*(factorial(k)^(c-1)))
                if Q ==0
                        println("yikes!")
                        return([1])
                else
                        Rone = div(r-1, Q)
                        Rtwo = r - Q*Rone
                        w = combination_unrank(Rone+1,n,k)
                        if n-k ==k
                                return([w, [1:k;]])
                        else
                                return([w, unrank_reg_combo(Rtwo, n-k, k)...])
                        end
                end
        end
end

function rank_reg_combo(v::Vector{Vector{Int}}, l::Int, r::Int)
        N= l*r
        s = 1
        for i in 1:(r-1)
                s = s + (combination_rank(v[i], N-l*(i-1), l)-1)*div(factorial(N-i*l), factorial(r-i)*factorial(l)^(r-i))
        end
        return(s)
end

function combo_class_next(v::Vector{Vector{Int}}, lambda::Vector{Int})
        n = sum(lambda)
        # Make dictionary assigning k to the number of k-cycles
        cc = counter(lambda)
        blocklengths = [k*cc[k] for k in keys(cc)]
        m = n
        outer_classes= [colex_bitstring(m -sum(blocklengths[1:(i-1)]), blocklengths[i]) for i in 1:length(blocklengths)]
        # Methods to consider - Split into blocks of k-cycles then do tuple iteration on a truncated version
        # Recurse "from the back"
end

function permutation_next(n, v::Vector{Int})
        k=n-1
        while v[k] > v[k+1]
                k = k-1
        end
        if k==0
                return(v)
        end
        j = n
        while v[k] >v[j]
                j = j-1
        end
        v[k],v[j] = v[j],v[k]
        r = n
        s = k+1
        while r >s
                v[r], v[s] = v[s], v[r]
                r = r-1
                s= s+1
        end
        return(v)
end

function colex_bitstring(n::Int,k::Int)#,v::Vector{Int})
        if k==0
                return([zeros(Int, n-k)])
        else
                if k < n
                        l1 = [vcat([0], C) for C in colex_bitstring(n-1, k)]
                        
                else
                        l1 = Vector{Int64}[]
                end
                l2 = [vcat([1], C) for C in colex_bitstring(n-1,k-1)]
                l = vcat(l1, l2)
                return(l)
        end
end

function conj_class_iter(n::Int, lambda::Vector{Int})
        #make dictionary assigns to each distinct value in lambda its multiplicity
        #l_mult = 

        #Chunk sizes - [l*mult(l) for each distinct l in lambda]

        #Construct default splitting and all perms within that splitting
        # We use combination_tuple_next on each chunk to get all the ways to split up the chunk and then go through the permutations from each split chunk
        # This will involve populating an array with a bunch of deepcopies.  We can initialize with appropriate size
        # Need to work out types of elements and number of elements

        #temporarily, rather than doings perms w/in a splitting, we'll just look at the partitions
#        N_combo_tuples = 
#        def_perms = 

        #Make list of chunkings (might need to recurse here)

        # Assign chunkings to each thread.  Within each chunking iterate of perm list and transform each to match the chunking

end

function combinations_list(n,k)
        V = colex_bitstring(n,k)
end

#function partition_type(v::Vector{Int})
#        
#end

function combination_rank(s::Vector{Int},n::Int, k::Int)
       return(binomial(n,k) - sum([binomial(n- s[i], k-i+1) for i in 1:k])) 
end

function combination_unrank(r::Int, n::Int, k::Int)
        if r > binomial(n,k)
                println("error r> n choose k")
                return([0])
        elseif n == k
                if r ==1
                        return [1:n;]
                else
                        println("error n=k, r not ")
                end
        elseif k == 1
                if r <= n
                        return [r]
                else
                        println("error")
                        return 0
                end
        else
                running_sum = binomial(n-1, k-1)
                i=1
                while running_sum <r
                        i += 1
                        running_sum += binomial(n-i, k-1)
                end
                running_sum = running_sum - binomial(n-i, k-1)
                out = vcat([i], [j+i for j in combination_unrank(r-running_sum, n-i, k-1)])
        end
        return(out)
end

function combination_unrank(k::Int,r::Int)
        v = zeros(Int64, k)
        for i in k:-1:1
                p = i-1
                while binomial(p,i) <= r
                        p = p+1
                end
                r = r - binomial(p-1, i)
                v[i] = p
        end
        return(v)
end

# Loop Plan: Store ways to partition into elements of each cycle-type
#  For first way, store all resulting permutations for those cycles in memory
#  Assign ways to partition to threads and loop over the permutation list, assign


# Computes the r-th partition of a set of size N into k blocks of size n
function regular_partition_unrank(r::Int, n::Int, k::Int)
        N = n*k
        total_parts = div(factorial(N), (factorial(n)^k)*factorial(k))
        fixed_first_block_num = div(total_parts, binomial(N-1, n-1))
        first_part_no = 1+ floor(r/fixed_first_block_num)
        new_r = r % fixed_first_block_num
end

struct regular_partition
        n::Int
        k::Int
        r::Int
end

#Base.length(P::regular_partition) = 
function Base.iterate(P::regular_partition)

end

#Method for special case of k 3-cycles in S_{3k}
#function conj_class_unrank(k::Int,r::Int)
#        n = 3*k
#        total_choices = div(factorial(n), factorial(k)*factorial(3)^k)
#        choices_for_one_cycle = n*(n-1)
#        one_cycles = [(1,i,j) for i != j, i,j in product(1:n)]
#        choices_per_cyc = total_choices/(n*(n-1))
#
#end

function ct_to_p(v::Vector{Vector{Int}})
        n = sum([length(w) for w in v])
        N = [1:n;]
        Q = [zeros(Int, length(w)) for w in v]
        for i in 1:length(v)
                Q[i] = N[v[i]]
                N= N[filter(x -> !(x in v[i]), eachindex(N))]
        end
        return(Q)
end
