include("search_organization.jl")
include("perm_utils.jl")
using Graphs

# Input should be the graphs found by make_minseps functions for genus g, along with the sets I_h for h < g
function make_I_g(g_graphlist, existing_graphs)
	I_g = []
	# for each graph in g_graphlist, want to iso test it againt all of existing_graphs, if it's not iso to any, it stays
	for graph in g_graphlist
		#println(string(graph))
		#flush(stdout)
		if graph == []
			println("hmm")
			flush(stdout)
			isDupe = 1
		else 
			isDupe = 0
			for e_graph in existing_graphs
                            if Graphs.Experimental.has_isomorph(graph, e_graph) ==1
					isDupe = 1
					break
				end
			end
			if isDupe == 0
				push!(I_g, graph)
			end
		end
	end
	return(I_g)
end

function main(E::Int)
    	graphmaketime = 0.0
    	graphisotime = 0.0
	I_g_list = [[] for i in 0:5]
       	I_gs = [Graphs.Graph([1;;])]
	I_g_list[1] = I_gs
	for g in 1:4
		#x = time()
		genus_g_duals= generate_minseps_genus(g)
		E_g_count = length(genus_g_duals)
		g_minseps = dual_list_to_minseps(genus_g_duals)
		println(string(g))
		#println(string(time()-x))
		flush(stdout)
		println("Size of E_g = ")
		println(string(E_g_count))
		flush(stdout)
       		#tt = time()
		graphsg = minseps_list_to_graphs(g_minseps, g)
       		#graphmaketime = graphmaketime + (time()-tt)
       		#tt = time()
		I_g = make_I_g(graphsg, I_gs)
		I_g_list[g+1] = I_g
		I_gs = vcat(I_gs, I_g)
       		#graphisotime = graphisotime + (time()-tt)
		println("Size of I_g = ")
		println(length(I_g))
		flush(stdout)
	end
	if E< 16
		g=5
		genus_g_duals = generate_minseps_genus(g,E)
		E_g_count = length(genus_g_duals)
		g_minseps = dual_list_to_minseps(genus_g_duals)
		println(string(g))
		flush(stdout)
		println("Size of E_g for $E edges")
		println(string(E_g_count))
		flush(stdout)
		graphsg = minseps_list_to_graphs(g_minseps, g)
		I_g = make_I_g(graphsg, I_gs)
		I_g_list[g+1] = I_g
		I_gs = vcat(I_gs, I_g)
		println("Size of I_g (for <= 15 Edges) = ")
		println(length(I_g))
		flush(stdout)
	elseif E ==16
		g =5
		total_minseps = []
		for ghat in [3,4,5]
			println("E=16, g=5, ghat = ")
			println(string(ghat))
			flush(stdout)
			push!(total_minseps, get_ghat_minseps_edges(g,ghat, E))
		end
		genus_g_duals = reduce(vcat, total_minseps)
		E_g_count = length(genus_g_duals)
		g_minseps = dual_list_to_minseps(genus_g_duals)
		println(string(g))
		flush(stdout)
		println("Size of E_g for 16 edges")
		println(string(E_g_count))
		flush(stdout)
		graphsg = minseps_list_to_graphs(g_minseps, g)
		I_g = make_I_g(graphsg, I_gs)
		println("Size of I_g for 16 edges")
		println(length(I_g))
		flush(stdout)
	elseif E<19
		g =5
		total_minseps = []
		for ghat in [4,5]
			println("E=$E, g=5, ghat = ")
			println(string(ghat))
			flush(stdout)
			push!(total_minseps, get_ghat_minseps_edges(g,ghat, E))
		end
		genus_g_duals = reduce(vcat, total_minseps)
		E_g_count = length(genus_g_duals)
		g_minseps = dual_list_to_minseps(genus_g_duals)
		println(string(g))
		flush(stdout)
		println("Size of E_g for $E edges")
		println(string(E_g_count))
		flush(stdout)
		graphsg = minseps_list_to_graphs(g_minseps, g)
		I_g = make_I_g(graphsg, I_gs)
		println("Size of I_g for $E edges")
		println(length(I_g))
		flush(stdout)
	else
		g =5
		println("E=$E, g=5, ghat = 5")
		#println(string(ghat))
		flush(stdout)
		total_minseps= get_ghat_minseps_edges(g,g, E)
		#genus_g_duals = reduce(vcat, total_minseps)
		E_g_count = length(total_minseps)
		g_minseps = dual_list_to_minseps(total_minseps)
		println(string(g))
		flush(stdout)
		println("Size of E_g for $E edges")
		println(string(E_g_count))
		flush(stdout)
		graphsg = minseps_list_to_graphs(g_minseps, g)
		println("graphslength 1 = ")
		println(string(length(graphsg)))
		flush(stdout)
		I_g = make_I_g(graphsg, I_gs)
		I_g_list[g+1] = I_g
		I_gs = vcat(I_gs, I_g)
		println("Size of I_g for $E edges")
		println(length(I_g))
		flush(stdout)
	end
#    println("total graphs make time = ")
#    println(graphmaketime)
#    println("total graphs iso time = ")
#    println(graphisotime)
end

main(13)
