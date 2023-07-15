include("august_make_minseps.jl")
include("check_transitive.jl")
include("graphmethods.jl")
using Graphs

# Input should be the graphs found by make_minseps functions for genus g, along with the sets I_h for h < g
function make_I_g(g_graphlist, existing_graphs)
	I_g = []
	# for each graph in g_graphlist, want to iso test it againt all of existing_graphs, if it's not iso to any, it stays
	for graph in g_graphlist
		if graph == []
			isDupe = 1
		else 
			isDupe = 0
			for e_graph in existing_graphs
				if are_isomorphic(graph, e_graph) ==1
                #if Graphs.Experimental.has_isomorph(Graphs.Graph(graph), e_graph) ==1
					isDupe = 1
					break
				end
			end
			if isDupe == 0
				#push!(I_g, Graphs.Graph(graph))
                push!(I_g, graph)
			end
		end
	end
	return(I_g)
end

function main(max_g::Int64)
    graphstime = 0.0
	I_g_list = [[] for i in 0:max_g]
	I_gs = [[1;;]]
	I_g_list[1] = I_gs
	for g in 1:max_g
		#x = time()
		genus_g_duals= generate_minseps_genus(g)
		g_minseps = dual_list_to_minseps(genus_g_duals)
		println(string(g))
		#println(string(time()-x))
		flush(stdout)
		println(string(length(g_minseps)))
		flush(stdout)
        tt = time()
		graphsg = minseps_list_to_graphs(g_minseps, g)
		I_g = make_I_g(graphsg, I_gs)
		I_g_list[g+1] = I_g
		I_gs = vcat(I_gs, I_g)
        graphstime= graphstime + time() - tt
		println(length(I_g))
		flush(stdout)
	end
    println("total graphs time = ")
    println(graphstime)
#	open("july_15_test_outputs.txt", "w") do file
#		for g in 0:max_g
#			write(file, string(g))
#			write(file, "\n")
#			for graph in I_g_list[g+1]
#				write(file, string(graph))
#				write(file, "\n")
#			end
#		end
#	end
end

main(1)
