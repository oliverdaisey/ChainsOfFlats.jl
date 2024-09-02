module ChainsOfFlats

using Oscar
numBitsConstant = 6 # number of bits to choose for generating random rationals

"""
    trivalent_tropical_linear_space(n::Int, k::Int, bergmanFan::Bool = false)

Given a number of variables `n` and number of linear forms `k`, returns a tropical linear space with exactly three maximal cones.
"""
function trivalent_tropical_linear_space(n::Int, k::Int, bergmanFan::Bool = false, numBits::Int = numBitsConstant)

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

    TropL = tropical_linear_space(I, nu)
    if (n_maximal_polyhedra(TropL) == 3)
        return TropL
    end
    # create the tropical linear space
    return trivalent_tropical_linear_space(n, k, bergmanFan, numBits)

end

end # module ChainsOfFlats
