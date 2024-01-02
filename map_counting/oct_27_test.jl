include("search_organization.jl")
include("perm_utils.jl")
using Graphs

function main()
	t = time()
	get_ghat_minseps_edges(5,5,14)
	tt = time()-t
	println(string(tt))
	println("done")
end

main()
