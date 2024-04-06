using PellsEquation
using Test

using .Iterators

@testset "PellsEquation.jl" begin

    @testset "pellsequation" begin
        pells13 = ((1, 0), (649, 180), (842401, 233640), (1093435849, 303264540))
        for (i, (x, y)) in zip(1:4, pellseqn(13))
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

        @test collect(take(pellseqn(991), 2)) == [(1, 0), (379516400906811930638014896080, 12055735790331359447442538767)]
        @test isempty(pellseqn(991, -1))

    end

    @testset "General Pell's Equation" begin
        @testset "x^2 - 157y^2 = ±12" begin
            ans = [
                (13, 1)
                (10663, 851)
                (579160, 46222)
                (483790960, 38610722)
                (26277068347, 2097138361)
                (21950079635497, 1751807067011)
                (1192216867392577, 95149264530709)
                (995897062658343427, 79481238398745359)
                (54092071464191542720, 4317017278925659678)
                (45184845607921619990920, 3606143265637668616178)
                (2454211373209617617356543, 195867390866986840719829)
                (2050081629081114757949687893, 163614326025765424385866679)
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

        @testset "x^2 - 5y^2 = ± 9" begin
            ans = [
                (3, 0)
                (27, 12)
                (483, 216)
                (8667, 3876)
                (155523, 69552)
                (2790747, 1248060)
                (50077923, 22395528)
                (898611867, 401871444)
                (16124935683, 7211290464)
                (289350230427, 129401356908)
                (5192179212003, 2322013133880)
                (93169875585627, 41666835052932)
            ]
            res = collect(take(pellseqn(5, 9), 12))
            @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
            @test res == ans

            ans = [
                (6, 3)
                (114, 51)
                (2046, 915)
                (36714, 16419)
                (658806, 294627)
                (11821794, 5286867)
                (212133486, 94868979)
                (3806580954, 1702354755)
                (68306323686, 30547516611)
                (1225707245394, 548152944243)
                (21994424093406, 9836205479763)
                (394673926435914, 176503545691491)
            ]
            res = collect(take(pellseqn(5, -9), 12))
            @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
            @test res == ans
        end

        @testset "x^2 - 20y^2 = ±80" begin
            ans = [
                (10, 1)
                (20, 4)
                (50, 11)
                (130, 29)
                (340, 76)
                (890, 199)
                (2330, 521)
                (6100, 1364)
                (15970, 3571)
                (41810, 9349)
                (109460, 24476)
                (286570, 64079)
            ]
            res = collect(take(pellseqn(20, 80), 12))
            @test res == ans

            ans = [
                (0, 2)
                (10, 3)
                (30, 7)
                (80, 18)
                (210, 47)
                (550, 123)
                (1440, 322)
                (3770, 843)
                (9870, 2207)
                (25840, 5778)
                (67650, 15127)
                (177110, 39603)
            ]
            res = collect(take(pellseqn(20, -80), 12))
            @test res == ans
        end
    end

    @testset "x^2 + 0 * y^2 = N^2" begin
        it = pellseqn(0, 3)
        @test eltype(it) == Tuple{BigInt,BigInt}
        @test isempty(it)

        it = pellseqn(0, 4)
        res = collect(Iterators.take(it, 5))
        @test eltype(it) == Tuple{BigInt,BigInt}
        @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
        @test res == [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4)]
    end

    @testset "x^2 + D * y^2 = 0" begin
        it = pellseqn(3, 0)
        @test eltype(it) == Tuple{BigInt,BigInt}
        @test collect(it) == [(0, 0)]

        it = pellseqn(4, 0)
        res = collect(Iterators.take(it, 5))
        @test eltype(it) == Tuple{BigInt,BigInt}
        @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
        @test res == [(0, 0), (2, 1), (4, 2), (6, 3), (8, 4)]
    end

    @testset "x^2 + D^2*y^2 = N" begin
        it = pellseqn(16, 10625)
        res = collect(Iterators.take(it, 5))
        @test eltype(it) == Tuple{BigInt,BigInt}
        @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
        @test isempty(it)
        @test res == [(105, 5), (225, 50), (321, 76), (1065, 265), (5313, 1328)]

        it = pellseqn(16, -9375)
        res = collect(Iterators.take(it, 6))
        @test eltype(it) == Tuple{BigInt,BigInt}
        @test all(xy isa Tuple{BigInt,BigInt} for xy in res)
        @test isempty(it)
        @test res == [(25, 25), (175, 50), (305, 80), (935, 235), (1561, 391), (4687, 1172)]
    end

    @testset "Continued fraction iterator" begin
        @testset "Continued fraction of 0" begin
            ans = Tuple{Int64,BigInt,BigInt}[(0, 0, 1)]
            @test collect(continued_fraction(0)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(0)) == ans
            @test collect(continued_fraction(0, 0)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(0, 0)) == ans
            @test collect(continued_fraction(0, 0, 1)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(0, 0, 1)) == ans
            @test_throws DomainError continued_fraction(0, 0, 0)
        end

        @testset "Continued fractions of integer" begin
            ans = Tuple{Int64,BigInt,BigInt}[(1, 1, 1)]
            @test collect(continued_fraction(1)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(1)) == ans
            @test collect(continued_fraction(0, 1)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(0, 1)) == ans
            @test collect(continued_fraction(0, 1, 1)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(0, 1, 1)) == ans
            @test collect(continued_fraction(1, 0, 1)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(1, 0, 1)) == ans

            ans = Tuple{Int64,BigInt,BigInt}[(16, 16, 1)]
            @test collect(continued_fraction(16^2)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(16^2)) == ans
            @test collect(continued_fraction(0, 16)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(0, 16)) == ans
            @test collect(continued_fraction(9, 13)) == ans
            @test collect(continued_fraction(36, 26, 2)) == ans

            ans = Tuple{Int64,BigInt,BigInt}[(-16, -16, 1)]
            @test collect(continued_fraction(16^2, 0, -1)) isa Vector{Tuple{Int64,BigInt,BigInt}}
            @test collect(continued_fraction(16^2, 0, -1)) == ans
            @test collect(continued_fraction(0, -16)) == ans
            @test collect(continued_fraction(1, -17)) == ans
            @test collect(continued_fraction(9, -51, 3)) == ans
        end
    end
end
