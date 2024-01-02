include("search_organization.jl")
include("perm_utils.jl")
using Graphs

function main()
	for i in 11:20
		t = time()
		println(string(i))
		Ecount = get_ghat_minseps_edges(5,5,i)
		println(string(Ecount))
		flush(stdout)
		tt = time()-t
		println(string(tt))
		println("done")
	end
end

main()
