module ChainsOfFlats

using Oscar
export trivalent_tropical_linear_space

numBitsConstant = 6 # number of bits to choose for generating random rationals
maxAttemptsConstant = 100 # maximum number of attempts to create a tropical linear space
verbose = true # print messages

"""
    trivalent_tropical_linear_space(n::Int, k::Int, bergmanFan::Bool = false, numBits::Int = numBitsConstant, maxAttempts::Int = maxAttemptsConstant)

Given positive integers `n` and `k`, this function attempts to create a trivalent tropical linear space in `n` variables with `k` linear forms.

The optional argument `bergmanFan` is a boolean that determines whether to perturb coefficients or not. The optional argument `numBits` is the number of bits to choose for generating random rationals. The optional argument `maxAttempts` is the maximum number of attempts to create a tropical linear space.
"""
function trivalent_tropical_linear_space(n::Int, k::Int, bergmanFan::Bool = false, numBits::Int = numBitsConstant, maxAttempts::Int = maxAttemptsConstant)

    # check if we have enough attempts left
    if (maxAttempts == 0)
        println("Ran out of attempts to create trivalent tropical linear space! Aborting.")
        return nothing
    end

    # the strategy is to create a linear ideal from binomials and one trinomial
    @assert (n > 0) & (k > 0) "The number of variables and linear forms must be positive."
    # @assert 2*k + 1 <= n "The number of variables must be at least two times the number of linear forms plus one."

    # create the polynomial ring
    Kt, t = rational_function_field(QQ, "t")
    nu = tropical_semiring_map(Kt, t, min)
    R, x = polynomial_ring(Kt, n)
    I = ideal([one(R)])

    while (saturation(I,ideal([prod(gens(R))])) == ideal([one(R)]))
        # create and push binomials
        linearForms = []
        for i in 1:(k-1)
            baseIndex = 2*i - 1
            randOne = rand(x)
            randTwo = rand(setdiff(x, [randOne]))
            f = rand(Int8)*randOne + rand(Int8)*randTwo
            push!(linearForms, f)
        end

        # create and push the trinomial
        randOne = rand(x)
        randTwo = rand(setdiff(x, [randOne]))
        randThree = rand(setdiff(x, [randOne, randTwo]))
        f = rand(Int8)*randOne + rand(Int8)*randTwo + rand(Int8)*randThree
        push!(linearForms, f)

        # for every form f, perturb the coefficients
        if !(bergmanFan)
            for i in 1:length(linearForms)
                linearForms[i] = map_coefficients(c -> c + rand_bits(QQ, numBits)*t^(rand_bits(ZZ, numBits)), linearForms[i])
            end
        end 
        
        # create the ideal
        I = ideal(linearForms)
    end

    if verbose println("Attempting to construct the tropical linear space with exactly three maximal cones.") end
    TropL = tropical_linear_space(I, nu)
    if (n_maximal_polyhedra(TropL) == 3)
        return TropL
    end

    if verbose println("Failed to create a tropical linear space with exactly three maximal cones. ($(maxAttempts-1) attempts left)") end

    return trivalent_tropical_linear_space(n, k, bergmanFan, numBits, maxAttempts-1)

end

end # module ChainsOfFlats
