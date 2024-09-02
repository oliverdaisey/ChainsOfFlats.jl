# Load the Test module
using Test
using ChainsOfFlats
using Oscar

# Run tests
@testset "ChainsOfFlats Tests" begin

    function trivalent_tropical_linear_space_test()
        println("Testing trivalent_tropical_linear_space")
        for pass in 1:8 # number of times to test
            n = rand(4:20)
            k = rand(1:(div(n, 2) - 1))
            println("Pass ", pass, " with n = ", n, " and k = ", k)
            T = ChainsOfFlats.trivalent_tropical_linear_space(n, k, false, 4)
            if n_maximal_polyhedra(polyhedral_complex(T)) != 3
                return false
            end
        end

        return true
    end


    @test trivalent_tropical_linear_space_test() == true

    
end