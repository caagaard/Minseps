include("perm_utils.jl")
#Evaluates the irreducible character of S_n corresponding to partition lambda for
#elements of cycle type rho.  Assumes that lambda has the for k^1,1^{n-k}.
function char_eval(n::Int, lambda::Vector{Int}, rho::Vector{Int})
        h = length(lambda)-1
        r = rho[1]
        charval = 0
        #edge cases:
        if r == n
                return((-1)^h)
        elseif length(lambda) ==1
                return(1)
        elseif length(lambda) ==n
                return((-1)^(r-1)*char_eval(n-r, ones(Int, n-r), rho[2:end]))
        end
        # contribution from leaving only horizontal strip
        if r< lambda[1]
                charval = charval + char_eval(n-r, vcat([lambda[1]-r], lambda[2:end]), rho[2:end])
        end
#        # contribution from leaving only vertical strip
        if r<length(lambda)
                charval = charval + (-1)^(r+1)*char_eval(n-r, vcat([lambda[1]], ones(Int, h-r)), rho[2:end])
        end
#
        return(charval)
end

function map_bound(n::Int, rho::Vector{Int})
        boundval = 0
        bigval = conj_class_size(rho)
        for i in 1:n
                lambda = vcat([i], ones(Int, n-i))
                boundval = boundval + bigval*char_eval(n, lambda, rho)/(n*char_eval(n, lambda, ones(Int, n)))
        end
        return(Int(round(boundval)))
end

function map_bound(n::Int, phi::Vector{Int}, theta::Vector{Int})
        boundval = 0
        bigval = conj_class_size(phi)*conj_class_size(theta)
        for i in 1:n
                lambda = vcat([i], ones(Int, n-i))
                boundval = boundval + bigval*char_eval(n,lambda, [n])*char_eval(n,lambda, phi)*char_eval(n,lambda,theta)/(factorial(big(n))*char_eval(n,lambda,ones(Int,n)))
        end
        return(Int(round(boundval)))
end
