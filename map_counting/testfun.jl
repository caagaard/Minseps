include("perm_utils.jl")

x = Perm([2,3,4,5,1])
println(string(conj_class_size(x)))
flush(stdout)
