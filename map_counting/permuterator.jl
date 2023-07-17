struct PermutationIterator{T}
    data::T
    length::Int
end

function has_repeats(state)
    for outer in 1:length(state)
        for inner in (outer+1):length(state)
            if state[outer] == state[inner]
                return true
            end
        end
    end
    return false
end

function increment!(state, max)
    state[end] +=1
    for i in length(state):-1:2
        if state[i] > max
            state[i] = 1
            state[i-1] += 1
        end
    end
end

function next_permutation!(state, max)
    while true
        increment!(state, max)
        if !has_repeats(state)
            break
        end
    end
end

function Base.iterate(p::PermutationIterator, state = ones(Int, p.length))
    next_permutation!(state, length(p.data))
    if state[1] > length(p.data)
        return nothing
    end
    [p.data[i] for i in state], state
end

Base.length(p::PermutationIterator) = Int(div(factorial(big(length(p.data))), factorial(big(length(p.data)-p.length))))
Base.eltype(p::PermutationIterator) = Vector{eltype(p.data)}

permutations(x, length) = PermutationIterator(x, length)
