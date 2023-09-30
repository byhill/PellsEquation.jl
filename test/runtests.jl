using PellsEquation
using Test

using .Iterators

@testset "PellsEquation.jl" begin

    @testset "pellsequation" begin
        pells13 = ((649, 180), (842401, 233640), (1093435849, 303264540))
        for (i, (x, y)) in zip(1:3, pellseqn(13))
            @test x isa BigInt
            @test y isa BigInt
            @test (x, y) == pells13[i]
        end

        negpells13 = ((18, 5), (23382, 6485), (30349818, 8417525))
        for (i, (x, y)) in zip(1:3, pellseqn(13, -1))
            @test x isa BigInt
            @test y isa BigInt
            @test (x, y) == negpells13[i]
        end

        @test first(pellseqn(991)) == (379516400906811930638014896080, 12055735790331359447442538767)
        @test isempty(pellseqn(991, -1))

    end

    @testset "General Pell's Equation" begin
        @testset "x^2 - 157y^2 = ±12" begin
            ans = [
                (13, 1),
                (10663, 851),
                (579160, 46222),
                (483790960, 38610722),
                (26277068347, 2097138361),
                (21950079635497, 1751807067011),
                (1192216867392577, 95149264530709),
                (995897062658343427, 79481238398745359),
                (54092071464191542720, 4317017278925659678),
                (45184845607921619990920, 3606143265637668616178),
                (2454211373209617617356543, 195867390866986840719829),
                (2050081629081114757949687893, 163614326025765424385866679),
            ]
            res = collect(Iterators.take(pellseqn(157, 12), 12))
            @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
            @test res == ans

            ans = [
                (50, 4)
                (2719, 217)
                (2271269, 181267)
                (123363799, 9845503)
                (103049745749, 8224265053)
                (5597138921710, 446700316396)
                (4675470012106610, 373143129538396)
                (253947789893540611, 20267240045357413)
                (212130749816239256561, 16929876922062299863)
                (11521865169662692139971, 919544947651210868827)
                (9624584245237121297322521, 768125445457745477545777)
                (522758544358818215189083630, 41720673799615848284192404)
            ]
            res = collect(Iterators.take(pellseqn(157, -12), 12))
            @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
            @test res == ans

        end
    end

end
