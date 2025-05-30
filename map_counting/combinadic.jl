using DataStructures

# Takes as input integers n and k, and k-subset of [n] = {1,2,...,n}
# Updatees v to be the next k-subset of [n] in lex order
function combination_next(n::Int,k::Int, v::Vector{Int})
        j=k
	# Checks if input is the last k-subset in lex order
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

# Given a vector of k-subsets of [N], [N-k], [N-2k],... computes the next such vector in lex order
function regular_combination_next!(v::Vector{Vector{Int}}, k::Int, N::Int)
        n = length(v)*k
        length_vec = [N:-k:(N-(length(v)-1)*k);]
        if v == [[1+N-n;length_vec[i]+2-k:length_vec[i];] for i in 1:length(length_vec)]
                return(v,1)
        end
        for i in length(v):-1:1
                if i == length(v)
                        firstval = 1+length_vec[i]-k
                else
                        firstval = 1+length_vec[i]-k*(1+length(v)-i)
                end
                if v[i] != [firstval; length_vec[i]+2-k:length_vec[i];]
                        v[i] = combination_next(length_vec[i],k, v[i])
                        for j in i+1:length(v)
                                v[j] = [v[j-1][1]: v[j-1][1]+k-1;]
                        end
                        return(v,0)
                end
        end
        return(v,0)
end

# Returns "building block for" the r-th set-partition of [n] where K gives the sizes of subsets in 
# the partition and mults[i] is the number of subsets of size K[i]
# Output should be fed into ct_to_p to obtain a permutation in the desired conjugacy class
function unrank_combo_partition(r::Int, n::Int, K::Vector{Int}, mults::Vector{Int})
        if length(K) != length(mults)
                println("K does not match mults")
                return([[0]])
        elseif length(K) ==1
                return(unrank_reg_combo(r,n,K[1]))
        end
        if n != sum([K[i]*mults[i] for i in 1:length(K)])
                println("error not a valid cycle type")
                return([[0]])
        elseif length(K)==1 && mults ==[1]
                return([[1:K[1];]])
        else
                qdenom =1
                for i in 2:length(K)
                        qdenom = qdenom * factorial(mults[i])*(factorial(K[i])^(mults[i]))
                end
                Q = div(factorial(n-(K[1]*mults[1])), qdenom)
                if Q ==0
                        println("yikes!")
                        return([[1]])
                else
                        # Here we're going to do a kind of converter.  We'll scale down by i-1 do a unrank_reg_combo, then scale back up by i-1
                        for i in 1:n-(K[1]*mults[1])
                                Rone = div(r-1, Q)
                                Rtwo = r-Q*Rone
                                flush(stdout)
                                w = unrank_reg_combo(Rone+1, n, K[1], mults[1])
                                return([w..., unrank_combo_partition(Rtwo, n-K[1]*mults[1], K[2:end], mults[2:end])...])
                        end
                end
        end
end

#This version unranks regular combinations, need to adapt to list of disjoint k-subsets
function unrank_reg_combo(r::Int, n::Int, k::Int)
        c = div(n,k)
        if c*k != n
                println("error k does not divide n")
                return([[0]])
        elseif n <=k
                return([[1:k;]])
        else
                Q = div(factorial(n-k), factorial(c-1)*(factorial(k)^(c-1)))
                if Q ==0
                        println("yikes!")
                        return([[1]])
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

function unrank_reg_combo(r::Int, n::Int, k::Int, l::Int)
        if l==0
                return(Vector{Int}[])
        elseif k*l > n
                println("error k*l > n")
                return([[0]])
        else
                for i in 1:(n+1-k*l)
                        Q = div(factorial(n+1-i-k), factorial(l-1)*factorial(k)^(l-1)*factorial(n+1-i-l*k))
                        if Q==0
                                println("yikes!")
                                return([[1]])
                        elseif r > Q*binomial(n-i, k-1)
                                r = r - Q*binomial(n-i,k-1)
                        else
                                Rone = div(r-1,Q)
                                Rtwo = r - Q*Rone
                                w = [i-1 for j in 1:k] + combination_unrank(Rone+1, n+1-i, k)
                                return([w, [[i-1 for j in 1:k] for m in 1:(l-1)]+unrank_reg_combo(Rtwo, n+1-i-k, k, l-1)...])
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

# input should be n (size of symmetric group) sigma a vector of vectors with entries ordered by decreasing length , and increasing first element within entries of the same length.  Within each 
# length block, the entries of the vectors should be integers in {1,... k} where k is the sum of the lengths of vectors in this block and all following blocks
# For now we're just going to iterate to the next cycle partition in the conjugacy class and handle the perm business later
function conj_class_next!(n::Int, sigma::Vector{Vector{Int}}, lambda::Vector{Int})
        chunk_track=0
        i = length(lambda)
        chunk_end = length(sigma)
        while i>0
                current_chunk = sigma[1+chunk_end - lambda[i] : chunk_end]
                t =  regular_combination_next!(current_chunk, length(current_chunk[1]), n - sum([length(sigma[i]) for i in 1:(chunk_end-lambda[i])]))[2]
                if t == 0
                        sigma[1+chunk_end - lambda[i] : chunk_end] = current_chunk
                        return(sigma, 0)
                else
                        sigma[1+chunk_end - lambda[i] : chunk_end] = [[1:length(current_chunk[1]);] for i in 1:lambda[i]]
                end
                chunk_end = chunk_end - lambda[i]
                i = i-1
        end
        return(sigma, 1)
end

function combination_rank(s::Vector{Int},n::Int, k::Int)
       return(binomial(n,k) - sum([binomial(n- s[i], k-i+1) for i in 1:k])) 
end

function combination_unrank(r::Int, n::Int, k::Int)
        running_sum = 1
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
                        return [0]
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

# Takes output of unrank_combo_partition and returns
# a "default" permutation corresponding to it
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
