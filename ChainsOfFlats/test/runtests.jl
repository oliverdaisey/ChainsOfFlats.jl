# Load the Test module
using Test
using ChainsOfFlats
using Oscar

# Run tests
@testset "ChainsOfFlats Tests" begin

    function trivalent_tropical_linear_space_test()
        for pass in 1:100
            println("Pass ", pass)
            n = rand(4:20)
            k = rand(1:(div(n, 2) - 1))
            T = ChainsOfFlats.trivalent_tropical_linear_space(n, k)
            if length(maximal_polyhedra(polyhedral_complex(T))) != 3
                return false
            end
        end

        return true
    end


    @test trivalent_tropical_linear_space_test() == true

    
end