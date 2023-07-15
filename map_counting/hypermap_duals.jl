include("perm_utils.jl")
using Oscar

function get_dual_map(psi::Perm, phi::Perm, n)
        phitemp = [n+phi[i] for i in Int8(1):n]
        theta = Perm(vcat([psi[i] for i in Int8(1):n],phitemp))
        alpha = Perm(vcat([i+n for i in Int8(1):n],[i for i in Int8(1):n]))
        sigma = theta^(-1)*alpha
        return([sigma,alpha])
end

function flipVert(vert, e)
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

function simplified_g(A = Matrix{Int4})
        (n,m) = size(A)
        for i in 1:m
                for j in i:n
                        while A[j,i] > 1
                                A[j,i] = A[j,i]-1
                                A[i,j] = A[j,i]
                                V = zeros(Int64, size(A)[1]+1)
                                V[j]=1
                                V[i] =1
                                A = hcat(A, V[1:end-1])
                                A = vcat(A,transpose(V))
                        end
                end
        end
        return(A)
end

function graph_from_embedding(nu, e, force_simple=0)
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
    if force_simple == 0
        return(Adj)
    else
        return(simplified_g(Adj))
    end
end 

