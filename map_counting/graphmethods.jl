include("permuterator.jl")
using IterTools
using LinearAlgebra
# Performance tests - array types for arrays of arrays and arrays of iterators
# Test at manageable size for each way to get perms

function makePerms(piG, v)
    if piG == [1]
        return(permutations([1:v;], v))
    else
        blocks = []
        firstblock = [1:piG[1]-1;]
        if length(firstblock) == 1
            push!(blocks, [firstblock])
        else
            push!(blocks, permutations(firstblock, length(firstblock)))
        end
        for i in [1:(length(piG)-1);]
            nextblock = [piG[i]:(piG[i+1]-1);]
            if length(nextblock) == 1
                push!(blocks, [nextblock])
            else
                push!(blocks, permutations(nextblock, length(nextblock)))
            end
        end
        endblock = [piG[end]:v;]
        if length(endblock) ==1
            push!(blocks, [endblock])
        else
            push!(blocks, permutations([piG[end]:v;], length([piG[end]:v;])))
        end
        # want product of set of perms of each block in blocks.
        tempperms = IterTools.product(blocks...)
    end
    speciallist = zip([1:length(tempperms);], tempperms)
    fullperms = Array{Array{Int64, 1},1}(undef, length(tempperms))
    for perm in speciallist
        fullperms[perm[1]] = vcat(perm[2]...)
    end
    return(fullperms)
end

function isConnected(Adj, v)
    if v==1
        return(true)
    elseif v==2
        if Adj[1,2] != 0
            return(true)
        else
            return(false)
        end
    else
        summation = Adj+Adj^0
        for i in 2:v-1
            summation = summation + Adj^i
        end
        if minimum(summation) ==0
            return(false)
        else
            return(true)
        end
    end
    #println("Yikes")
    return(false)
end

function getblocks(v)
    blocks = Int64[]
    for i in [2:length(v);]
        if v[i]!= v[i-1]
            push!(blocks, i)
        end
    end
    if length(blocks) == 0
        push!(blocks,1)
    end
    return(blocks)
end

function default_makeCanon(A)
    ranking= (vec(sum(A, dims=1))+diag(A)).*200 + diag(A).*20+count.(>(0), eachcol(A))
    M = A[sortperm(ranking), sortperm(ranking)]
    sortedrank = sort(ranking)
    return(M, sortedrank, getblocks(sortedrank))
end

function makeAdj(mat)
    rows = split(mat[2:end-1], "; ")
    Adjmat = parse.(Int64, split(rows[1], " "))
    for row in rows[2:end]
        Adjmat = hcat(Adjmat, parse.(Int64, split(row, " ")))
    end
    return(Adjmat)
end

# Returns whether two graphs are isomorphic
function are_isomorphic(graph1, graph2)
	G1 = default_makeCanon(graph1)
	G2 = default_makeCanon(graph2)
	if G1[2] == G2[2]
		v = size(graph1,1)
		permList = makePerms(G1[3], v)
		if length(permList) == 1
			if G1[1] == G2[1]
				return(1)
			else
				return(0)
			end
		else
			permsOfTarget = [G1[1][sortperm(permVec),sortperm(permVec)] for permVec in permList]
			for perm in permsOfTarget
				if perm == G2[1]
					return(1)
				end
			end
		end
	end
	return(0)
end

function getIsoClasses(graphlist::Vector{Matrix{Int64}})
	iso_classes = []
	just_graphs = Vector{Matrix{Int64}}()
	just_vecs = []
#	v = size(graphlist[1], 1)
	for A in graphlist
		v= size(A,1)
		G = default_makeCanon(A)
		(D, U) = (G[1] in just_graphs, G[2] in just_vecs)
		if  D
			ghostvar =1;
		else
			if !U
				push!(just_graphs, G[1])
				push!(just_vecs, G[2])
				push!(iso_classes, G)
			else
				permList = makePerms(G[3], v)
				permsOfTarget = []
				for permVec in permList
					push!(permsOfTarget, G[1][sortperm(permVec), sortperm(permVec)])
				end
				isDupe = 0
				for testgraph in iso_classes
					if G[2] == testgraph[2]
						# Do stuff
						for perm in permsOfTarget
							if perm == testgraph[1]
								isDupe =1
								break
							end
						end
					end
					if isDupe ==1
						break
					end
				end 
				if isDupe == 0
					push!(just_graphs, G[1])
					push!(just_vecs, G[2])
					push!(iso_classes, G)
				end
			end
		end
	end
	return(just_graphs)
end

#function getIsoClasses(graphlist::Array{Any, 2}, e)
#	just_graphs = []
#	just_vecs = []
#	v = size(graphlist[1,1], 1)
#	# First we initialize iso_classes
#	A = graphlist[1,:]
#	G = default_makeCanon(A[1])
#	push!(just_graphs, G[1])
#	push!(just_vecs, G[2])
#	iso_classes = [G A[2]]	
#	# Now we loop over the rest of the graphs
#	for i in 2:size(graphlist, 2)
#		A = graphlist[i,:]
#		#println(A)
#		#println(A[1])
#		G = default_makeCanon(A[1])
#		(D, U) = (G[1] in just_graphs, G[2] in just_vecs)
#		if  D
#			ghostvar =1;
#		else
#			if !U
#				push!(just_graphs, G[1])
#				push!(just_vecs, G[2])
#				iso_classes = vcat(iso_classes, [G A[2]])
#			else
#				permList = makePerms(G[3], v)
#				permsOfTarget = []
#				for permVec in permList
#					push!(permsOfTarget, G[1][sortperm(permVec), sortperm(permVec)])
#				end
#				isDupe = 0
#				for testgraph in iso_classes[:,1]
#					if G[2] == testgraph[2]
#						# Do stuff
#						for perm in permsOfTarget
#							if perm == testgraph[1]
#								isDupe =1
#								break
#							end
#						end
#					end
#					if isDupe ==1
#						break
#					end
#				end 
#				if isDupe == 0
#					push!(just_graphs, G[1])
#					push!(just_vecs, G[2])
#					iso_classes = vcat(iso_classes, [G A[2]])
#				end
#			end
#		end
#	end
#	#print("0")
#	for x in iso_classes
#	#	println(string(x))
#	end
#	graph_system_classes = [[grph[1][1], grph[2]] for grph in eachrow(iso_classes)]
#	if length(graph_system_classes) != length(just_graphs)
#		println("Oh No, new method is messed up")
#	end
#	return(graph_system_classes)
#end

