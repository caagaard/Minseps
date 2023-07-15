using AbstractAlgebra

function find_partitions(e::Int, n::Int)
	if n == 1
		return([e])
	elseif e == n
		return([[1 for i in 1:e]])
	else
		partitions = []
		for i in (e-n+1):-1:div((e+n-1), n)
			top_block = [i]
			subpartitions = find_partitions(e-i, n-1)
			tempparts = [vcat(top_block, subpartition) for subpartition in subpartitions]
			partitions = vcat(partitions, tempparts)
		end
	end
	return(partitions)
end

function lexVal(chunk)
	lexSum =0
	for i in 1:length(chunk)
		lexSum += chunk[i]*(10^(length(chunk)-i))
	end
	return(lexSum)
end	

function isOrderedPart(parttuple)
	M = [lexVal(chunk) for chunk in parttuple]
	return(issorted(M, rev=true))
end

function get_partitions(e, big_blocks)
	# Need to figure out ways to partition each side
	component_one= find_partitions(e,big_blocks[1])
	component_two = find_partitions(e, big_blocks[2])
	if big_blocks[1] != big_blocks[2]
		return(Iterators.product(component_one, component_two))
	else
		temp_parts = [[pair, isOrderedPart(pair)] for  pair in Iterators.product(component_one, component_two)]
		tobezipped_1 = [pair[1][1] for pair in temp_parts if pair[2] == 1]
		tobezipped_2 = [pair[1][2] for pair in temp_parts if pair[2] == 1]
		return(Iterators.zip(tobezipped_1,tobezipped_2))
	end
	# Then need to combine somehow
	# Then worry about duplicates
end	

function build_theta(e, stru)
	perm = Int8[]
	starter = Int8(0)
	for block in stru
		for x in block
			starter = starter +x
			perm = vcat(perm, [starter; (starter-x+1):(starter -1);])
		end
	end
	return(perm)
end

function isFastTheta(poss_struct,e)
	if poss_struct[1] == [e]
		return(1)
	else
		overlap_compare_1 = countmap(poss_struct[1])
		overlap_compare_2 = countmaps(poss_struct[2])
		for (i, j) in Iterators.product(keys(overlap_compare_1), keys(overlap_compare_2))
			if i*overlap_compare_1[i]+j*overlap_compare_2[j] > e
				return(1)
			end
		end
		return(2)
	end
end

function makethetas(e, g, v)
	nsums = v+2+2*g-e
	boundarydiv = [[i,nsums-i] for i in 1:div(nsums,2)]
	thetas = []
	fast_thetas = []
	slow_thetas = []
	for outline in boundarydiv
#		println(outline)
		poss_structs = get_partitions(e, outline)
		for poss_struct in poss_structs
#			is_fast_poss_struct = isFastTheta(poss_struct, e)
			theta = build_theta(e, poss_struct)
			push!(thetas, Perm(theta))
		end
	end
	return(thetas)
end
