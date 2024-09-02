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
    @assert 2*k + 1 <= n "The number of variables must be at least two times the number of linear forms plus one."

    # create the polynomial ring
    Kt, t = rational_function_field(QQ, "t")
    nu = tropical_semiring_map(Kt, t, min)
    R, (x...) = polynomial_ring(Kt, ["x$i" for i in 1:n])
    x=x[1]

    # create and push binomials
    linearForms = []
    for i in 1:(k-1)
        baseIndex = 2*i - 1
        f = x[baseIndex] + x[baseIndex + 1]
        push!(linearForms, f)
    end

    # create and push the trinomial
    f = x[2*k - 1] + x[2 * k] + x[2*k + 1]
    push!(linearForms, f)

    # for every form f, perturb the coefficients
    if !(bergmanFan)
        for i in 1:length(linearForms)
            linearForms[i] = map_coefficients(c -> c + rand_bits(QQ, numBits)*t^(rand_bits(ZZ, numBits)), linearForms[i])
        end
    end 
    
    # create the ideal
    I = ideal(linearForms)

    # create the tropical linear space
    return tropical_linear_space(I, nu)

end

end # module ChainsOfFlats
