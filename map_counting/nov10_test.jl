include("search_organization.jl")
include("perm_utils.jl")
using Graphs

function main()
	Ecount = get_ghat_minseps_edges(5,5,15)
	println(string(Ecount))
	flush(stdout)
end

main()
