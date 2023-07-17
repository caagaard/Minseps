using GAP
using Oscar

function cycle_type(phi)
        v= [length(x) for x in cycles(phi)]
        return(sort(v))
end

# Takes a pair of permutations of type Perm{IntN} for some N and returns whether they generate a transitive permutation group
function is_trans(permpair::Vector{Perm{Int8}})
        phi = permpair[1]
        psi = permpair[2]
        n = sum(cycle_type(phi))
        phi = GAP.Globals.PermList(GapObj([Int64(phi[i]) for i in 1:n]))
        psi = GAP.Globals.PermList(GapObj([Int64(psi[i]) for i in 1:n]))
        H = GAP.Globals.Group(phi,psi)
        return(GAP.Globals.IsTransitive(H))
end

# Takes an integer 1 and permutations phi, psi.
# Returns [n,m], where n is the length of the obit of i under phi, and m is the length of the orbit of i under psi
function point_type(i, phi,psi)
        return([[length(cycles(phi)[j]) for j in 1:length(cycles(phi)) if i in cycles(phi)[j]][1], [length(cycles(psi)[j]) for j in 1:length(cycles(psi)) if i in cycles(psi)[j]][1]])
end

# Takes a pair of integers x and y, and permutations phi and psi.  Returns true if there is a color swapping isomorphism mapping x to y
function is_color_swap_iso(x,y,phi,psi)
        E = sum(cycle_type(phi))
        f_array = [0 for i in 1:E]
        # Fill in f_array
        f_array[x] = y
        counter = 0
        while minimum(f_array)== 0 && counter <= E
                for i in 1:E
                        if f_array[i] != 0
                                f_array[phi[i]] = psi[f_array[i]]
                                f_array[psi[i]] = phi[f_array[i]]
                        end
                end
                counter = counter+1
        end
        if sort(f_array) != [1:E;]
                return(0)
        end
        f = Perm(f_array)
                if f*phi == psi*f && f*psi == phi*f
                        return(1)
                end
        return(0)
end

# Determines if a hypermap is `self-color-dual'
# Input should be a pair of permutations
function is_self_color_dual(hypermap)
        phi = hypermap[1]
        psi = hypermap[2]
        E = sum(cycle_type(phi))
        # first easy check is if cycle structures match
        if cycle_type(phi) != cycle_type(psi)
                return(0)
        end
        # Now we check for automorphism:
        # Possible destinations for 1: set n = length of 1s cycle in phi, m = lengh of 1s cycle in psi
        # Destinations are all other elements with type (m,n)
        one_targets = [y for y in 2:E if point_type(y,phi,psi) == point_type(1, psi, phi)]
        for y in one_targets
                if(is_color_swap_iso(1,y,phi,psi)) == 1
                        return(1)
                end
        end
        return(0)
end
