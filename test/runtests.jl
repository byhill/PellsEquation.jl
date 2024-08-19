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

        T = NTuple{3,BigInt}
        f(itr) = collect(take(itr, 10))

        @testset "Invalid continued fraction inputs" begin
            @test_throws DomainError continued_fraction(0, 0, 0)
            @test_throws DomainError continued_fraction(7890789132, -13458678056781065478608, 0)
            @test_throws DomainError continued_fraction(-1)
            @test_throws DomainError continued_fraction(-1, 0, 0)
            @test_throws DomainError continued_fraction(typemin(Int128), 3213789, 19876)
        end

        @testset "Continued fractions of integers" begin
            @testset "Continued fraction of 0" begin
                ans = T[(0, 0, 1)]
                @test collect(continued_fraction(0)) isa Vector{T}
                @test collect(continued_fraction(0)) == ans
                @test collect(continued_fraction(0, 0)) isa Vector{T}
                @test collect(continued_fraction(0, 0)) == ans
                @test collect(continued_fraction(0, 0, 1)) isa Vector{T}
                @test collect(continued_fraction(0, 0, 1)) == ans
                @test collect(continued_fraction(9, -3, 3)) isa Vector{T}
                @test collect(continued_fraction(9, -3, 3)) == ans
            end

            @testset "Continued fraction of 1" begin
                ans = T[(1, 1, 1)]
                @test collect(continued_fraction(1)) isa Vector{T}
                @test collect(continued_fraction(1)) == ans
                @test collect(continued_fraction(1, 0)) isa Vector{T}
                @test collect(continued_fraction(1, 0)) == ans
                @test collect(continued_fraction(0, 1)) isa Vector{T}
                @test collect(continued_fraction(0, 1)) == ans
                @test collect(continued_fraction(0, 1, 1)) isa Vector{T}
                @test collect(continued_fraction(0, 1, 1)) == ans
                @test collect(continued_fraction(1, 0, 1)) isa Vector{T}
                @test collect(continued_fraction(1, 0, 1)) == ans
            end

            @testset "Continued fraction of -1" begin
                ans = T[(-1, -1, 1)]
                @test collect(continued_fraction(1, 0, -1)) isa Vector{T}
                @test collect(continued_fraction(1, 0, -1)) == ans
                @test collect(continued_fraction(0, -1)) isa Vector{T}
                @test collect(continued_fraction(0, -1)) == ans
                @test collect(continued_fraction(1, -2)) isa Vector{T}
                @test collect(continued_fraction(1, -2)) == ans
                @test collect(continued_fraction(4, -4, 2)) isa Vector{T}
                @test collect(continued_fraction(4, -4, 2)) == ans
            end

            @testset "Continued fraction of 16" begin
                ans = T[(16, 16, 1)]
                @test collect(continued_fraction(16^2)) isa Vector{T}
                @test collect(continued_fraction(16^2)) == ans
                @test collect(continued_fraction(0, 16)) == ans
                @test collect(continued_fraction(9, 13)) == ans
                @test collect(continued_fraction(36, 26, 2)) == ans
            end

            @testset "Continued fraction of -16" begin
                ans = T[(-16, -16, 1)]
                @test collect(continued_fraction(16^2, 0, -1)) isa Vector{T}
                @test collect(continued_fraction(16^2, 0, -1)) == ans
                @test collect(continued_fraction(0, -16)) == ans
                @test collect(continued_fraction(1, -17)) == ans
                @test collect(continued_fraction(9, -51, 3)) == ans
            end
        end

        @testset "Continued fractions of rational numbers" begin
            @testset "Continued fraction of 20 / 3" begin
                ans = T[(6, 6, 1), (1, 7, 1), (2, 20, 3)]
                @test collect(continued_fraction(0, 20, 3)) isa Vector{T}
                @test collect(continued_fraction(0, 20, 3)) == ans
                @test collect(continued_fraction(0, -20, -3)) == ans
                @test collect(continued_fraction(1, -21, -3)) == ans
            end

            @testset "Continued fraction of -20 / 3" begin
                ans = T[(-7, -7, 1), (3, -20, 3)]
                @test collect(continued_fraction(0, -20, 3)) isa Vector{T}
                @test collect(continued_fraction(0, -20, 3)) == ans
                @test collect(continued_fraction(0, 20, -3)) == ans
                @test collect(continued_fraction(4, 18, -3)) == ans
            end

            @testset "Continued fraction of 234987 / 32482" begin
                ans = T[(7, 7, 1), (4, 29, 4), (3, 94, 13), (1, 123, 17), (3, 463, 64), (253, 117262, 16209), (2, 234987, 32482)]
                @test collect(continued_fraction(0, 234987, 32482)) isa Vector{T}
                @test collect(continued_fraction(0, 234987, 32482)) == ans
                @test collect(continued_fraction(
                    153238977834008140323061335508155135750477250074142538975866730199164956756869883813272310084,
                    29135547569801845033916408910314288308671267275763353765620852976,
                    4027375370391994156245948082929659444081323974680946008264032844
                )) == ans
            end

            @testset "Continued fraction of -234987 / 32482" begin
                ans = [(-8, -8, 1), (1, -7, 1), (3, -29, 4), (3, -94, 13), (1, -123, 17), (3, -463, 64), (253, -117262, 16209), (2, -234987, 32482)]
                @test collect(continued_fraction(0, -234987, 32482)) isa Vector{T}
                @test collect(continued_fraction(0, -234987, 32482)) == ans
                @test collect(continued_fraction(
                    54802076066125638907663991125841439074608825953610000000000000000,
                    -231953341023216900028432549687813208662698133100000000,
                    32062660585973400003897841555819519302932860000000000
                )) == ans
            end
        end

        @testset "Continued fractions of quadratic integers" begin
            @testset "Continued fraction of sqrt(3)" begin
                ans = T[(1, 1, 1), (1, 2, 1), (2, 5, 3), (1, 7, 4), (2, 19, 11), (1, 26, 15), (2, 71, 41), (1, 97, 56), (2, 265, 153), (1, 362, 209)]
                @test f(continued_fraction(3)) isa Vector{T}
                @test f(continued_fraction(3)) == ans
                @test f(continued_fraction(3, 0)) == ans
                @test f(continued_fraction(3, 0, 1)) == ans
            end

            @testset "Continued fraction of -sqrt(3)" begin
                ans = [(-2, -2, 1), (3, -5, 3), (1, -7, 4), (2, -19, 11), (1, -26, 15), (2, -71, 41), (1, -97, 56), (2, -265, 153), (1, -362, 209), (2, -989, 571)]
                @test f(continued_fraction(3, 0, -1)) isa Vector{T}
                @test f(continued_fraction(3, 0, -1)) isa Vector{T}
                @test f(continued_fraction(116964995676226639042340992193288944844, 0, -12345678901234567890)) isa Vector{T}
            end

            @testset "Continued fraction of sqrt(3) - 1" begin
                ans = [(0, 0, 1), (1, 1, 1), (2, 2, 3), (1, 3, 4), (2, 8, 11), (1, 11, 15), (2, 30, 41), (1, 41, 56), (2, 112, 153), (1, 153, 209)]
                @test f(continued_fraction(3, -1, 1)) isa Vector{T}
                @test f(continued_fraction(3, -1, 1)) == ans
                @test f(continued_fraction(12, -2, 2)) == ans
            end

            @testset "Continued fraction of -sqrt(3) + 1" begin
                ans = T[(-1, -1, 1), (3, -2, 3), (1, -3, 4), (2, -8, 11), (1, -11, 15), (2, -30, 41), (1, -41, 56), (2, -112, 153), (1, -153, 209), (2, -418, 571)]
                @test f(continued_fraction(3, -1, -1)) isa Vector{T}
                @test f(continued_fraction(3, -1, -1)) == ans
                @test f(continued_fraction(12, -2, -2)) == ans
            end

            @testset "Continued fraction of sqrt(3) - 2" begin
                ans = T[(-1, -1, 1), (1, 0, 1), (2, -1, 3), (1, -1, 4), (2, -3, 11), (1, -4, 15), (2, -11, 41), (1, -15, 56), (2, -41, 153), (1, -56, 209)]
                @test f(continued_fraction(3, -2)) isa Vector{T}
                @test f(continued_fraction(3, -2)) == ans
                @test f(continued_fraction(48, -8, 4)) == ans
            end

            @testset "Continued fraction of -sqrt(3) + 2" begin
                ans = T[(0, 0, 1), (3, 1, 3), (1, 1, 4), (2, 3, 11), (1, 4, 15), (2, 11, 41), (1, 15, 56), (2, 41, 153), (1, 56, 209), (2, 153, 571)]
                @test f(continued_fraction(3, -2, -1)) isa Vector{T}
                @test f(continued_fraction(3, -2, -1)) == ans
                @test f(continued_fraction(27, -6, -3)) == ans
            end
        end

        @testset "Continued fractions of quadratic irrationals" begin

            @testset "Continued fraction of (sqrt(55) - 92031) / 26" begin
                ans = T[(-3540, -3540, 1), (1, -3539, 1), (1, -7079, 2), (1, -10618, 3), (2, -28315, 8), (2, -67248, 19), (14, -969787, 274), (2, -2006822, 567), (2, -4983431, 1408), (2, -11973684, 3383)]
                @test f(continued_fraction(55, -92031, 26)) isa Vector{T}
                @test f(continued_fraction(55, -92031, 26)) == ans
            end

            @testset "Continued fraction of (sqrt(1237897) - 12309821389031231) / -4238973421237" begin
                ans = T[(0, 0, 1), (344357020891907793, 1, 344357020891907793), (1, 1, 344357020891907794), (7, 8, 2754856167135262351), (1, 9, 3099213188027170145), (12, 116, 39945414423461304091), (2, 241, 82990042034949778327), (1, 357, 122935456458411082418), (2, 955, 328860954951771943163), (1, 1312, 451796411410183025581)]
                @test f(continued_fraction(1237897, -12309821389031231, -4238973421237897897893123123123543)) isa Vector{T}
                @test f(continued_fraction(1237897, -12309821389031231, -4238973421237897897893123123123543)) == ans
            end

            @testset "Continued fraction of (sqrt(8963213689631268963232167896786742675) + 3) / 4" begin
                ans = T[(748465667617396002, 748465667617396002, 1), (1, 748465667617396003, 1), (2, 2245397002852188008, 3), (2, 5239259673321772019, 7), (4, 23202435696139276084, 31), (1, 28441695369461048103, 38), (4, 136969217173983468496, 183), (3, 439349346891411453591, 587), (1, 576318564065394922087, 770), (70, 40781648831469055999681, 54487)]
                @test f(continued_fraction(8963213689631268963232167896786742675, 3, 4)) isa Vector{T}
                @test f(continued_fraction(8963213689631268963232167896786742675, 3, 4)) == ans
            end
        end

        @testset "Continued fraction using different input types" begin

            @testset "Continued fraction of (sqrt(3) + 4) / 5 using different input types" begin
                ans = T[(1, 1, 1), (6, 7, 6), (1, 8, 7), (4, 39, 34), (1, 47, 41), (7, 368, 321), (1, 415, 362), (4, 2028, 1769), (1, 2443, 2131), (7, 19129, 16686)]
                for U in [Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32, UInt64, UInt128, BigInt]
                    @testset "Continued fraction of (sqrt(3) + 4) / 5 with input type $U" begin
                        @test f(continued_fraction(U(3), U(4), U(5))) isa Vector{T}
                        @test f(continued_fraction(U(3), U(4), U(5))) == ans
                    end
                end
            end

            @testset "Continued fraction of (sqrt(3) - 4) / 5 using different input types" begin
                ans = T[(-1, -1, 1), (1, 0, 1), (1, -1, 2), (4, -4, 9), (1, -5, 11), (7, -39, 86), (1, -44, 97), (4, -215, 474), (1, -259, 571), (7, -2028, 4471)]
                for U in [Int8, Int16, Int32, Int64, Int128, BigInt]
                    @testset "Continued fraction of (sqrt(3) - 4) / 5 with input type $U" begin
                        @test f(continued_fraction(U(3), U(-4), U(5))) isa Vector{T}
                        @test f(continued_fraction(U(3), U(-4), U(5))) == ans
                    end
                end
            end
        end
    end
end
