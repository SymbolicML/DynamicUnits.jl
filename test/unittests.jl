using DynamicQuantities
using DynamicQuantities: NoDims, AbstractSymbolicDimensions
using DynamicQuantities: DEFAULT_QUANTITY_TYPE, DEFAULT_DIM_BASE_TYPE, DEFAULT_DIM_TYPE, DEFAULT_VALUE_TYPE
using DynamicQuantities: array_type, value_type, dim_type, quantity_type
using DynamicQuantities: GenericQuantity, with_type_parameters, constructorof
using DynamicQuantities: promote_quantity_on_quantity, promote_quantity_on_value
using DynamicQuantities: UNIT_VALUES, UNIT_MAPPING, UNIT_SYMBOLS, ALL_MAPPING, ALL_SYMBOLS, ALL_VALUES
using DynamicQuantities.SymbolicUnits: SYMBOLIC_UNIT_VALUES
using DynamicQuantities.AffineUnits: AFFINE_UNIT_SYMBOLS, AFFINE_UNIT_MAPPING, AFFINE_UNIT_VALUES
using DynamicQuantities: map_dimensions
using DynamicQuantities: _register_unit, _register_affine_unit
using Ratios: SimpleRatio
using SaferIntegers: SafeInt16
using StaticArrays: SArray, MArray
using LinearAlgebra: norm
using Test

function record_show(s, f=show)
    io = IOBuffer()
    f(io, s)
    return String(take!(io))
end
function unsafe_isapprox(x, y; kwargs...)
    return isapprox(ustrip(x), ustrip(y); kwargs...) && dimension(x) == dimension(y)
end

# TODO: This is a bit hacky but is required to avoid ambiguities
Base.round(::Type{T}, x::SimpleRatio) where {T} = round(T, x.num // x.den)

@testset "Basic utilities" begin

    for Q in [Quantity, GenericQuantity, RealQuantity], T in [DEFAULT_VALUE_TYPE, Float16, Float32, Float64], R in [DEFAULT_DIM_BASE_TYPE, Rational{Int16}, Rational{Int32}, SimpleRatio{Int}, SimpleRatio{SafeInt16}]
        D = Dimensions{R}
        x = Q(T(0.2), D, length=1, mass=2.5)

        @test typeof(x).parameters[1] == T
        @test typeof(x).parameters[2] == D
        @test ulength(x) == R(1 // 1)
        @test umass(x) == R(5 // 2)
        @test ustrip(x) ≈ T(0.2)
        @test dimension(x) == Dimensions(R, length=1, mass=5 // 2)
        if R == DEFAULT_DIM_TYPE
            @test x == Q(T(0.2), length=1, mass=2.5)
            @test dimension(x) == Dimensions(length=1, mass=5 // 2)
        end

        y = x^2

        @test typeof(x).parameters[1] == T
        @test typeof(x).parameters[2] == D
        @test ulength(y) == R(2 // 1)
        @test umass(y) == (5 // 1)
        @test ustrip(y) ≈ T(0.04)

        y = x + x

        @test ulength(y) == R(1 // 1)
        @test umass(y) == R(5 // 2)
        @test ustrip(y) ≈ T(0.4)

        if R <: Rational
            @test string(x) == "0.2 m kg⁵ᐟ²"
            @test string(inv(x)) == "5.0 m⁻¹ kg⁻⁵ᐟ²"
        end

        @test_throws DimensionError x^2 + x

        # Check output of error:
        try
            x^2 + x
            @test false
        catch e
            r"DimensionError: .* and .* have incompatible dimensions."
            io = IOBuffer()
            showerror(io, e)
            msg = String(take!(io))
            @test occursin("DimensionError", msg)
            @test occursin("incompatible dimensions", msg)
        end

        y = inv(x)

        @test ulength(y) == R(-1 // 1)
        @test umass(y) == R(-5 // 2)
        @test ustrip(y) ≈ R(5)

        y = x - x

        @test iszero(x) == false
        @test iszero(y) == true
        @test iszero(y.dimensions) == false

        y = -x

        @test ustrip(y) == -ustrip(x)
        @test dimension(y) == dimension(x)

        y = x / x

        @test iszero(x.dimensions) == false
        @test iszero(y.dimensions) == true

        y = Q(T(2 // 10), D, length=1, mass=5 // 2)

        @test y ≈ x

        y = Q(T(2 // 10), D, length=1, mass=6 // 2)

        @test !(unsafe_isapprox(y, x))

        y = x * Inf32

        @test typeof(y).parameters[1] == promote_type(T, Float32)
        @test typeof(y).parameters[2] == D
        @test isfinite(x)
        @test !isfinite(y)

        y = x^2.1

        @test typeof(y).parameters[1] == T  # Should not promote! Expect 2.1 to be converted to 21//10
        @test typeof(y).parameters[2] == D
        @test ulength(y) == R(1 * (21 // 10))
        @test umass(y) == R((5 // 2) * (21 // 10))
        @test utime(y) == R(0)
        @test ucurrent(y) == R(0)
        @test utemperature(y) == R(0)
        @test uluminosity(y) == R(0)
        @test uamount(y) == R(0)
        @test ustrip(y) ≈ T(0.2^2.1)

        dimensionless = Q(one(T), D)
        y = T(2) + dimensionless
        @test ustrip(y) == T(3)
        @test dimension(y) == Dimensions(R)
        @test typeof(y) == Q{T,D}

        y = T(2) - dimensionless
        @test ustrip(y) == T(1)
        @test dimension(y) == Dimensions(R)
        @test typeof(y) == Q{T,D}

        y = dimensionless + T(2)
        @test ustrip(y) == T(3)
        y = dimensionless - T(2)
        @test ustrip(y) == T(-1)

        @test_throws DimensionError Q(one(T), D,  length=1) + 1.0
        @test_throws DimensionError Q(one(T), D, length=1) - 1.0
        @test_throws DimensionError 1.0 + Q(one(T), D, length=1)
        @test_throws DimensionError 1.0 - Q(one(T), D, length=1)
    end

    x = Quantity(-1.2, length=2 // 5)

    @test typemax(x) == Quantity(typemax(-1.2), length=2 // 5)
    @test typemax(typeof(x)) == Quantity(typemax(typeof(-1.2)))

    @test abs(x) == Quantity(1.2, length=2 // 5)
    @test abs(x) == abs(Quantity(1.2, length=2 // 5))
    @test abs2(x) == Quantity(abs2(-1.2), length=4 // 5)

    @test copy(x) == x

    @test iszero(x) == false
    @test iszero(x * 0) == true
    @test isfinite(x) == true
    @test isfinite(x * Inf) == false
    @test isfinite(x * NaN) == false
    @test isinf(x * Inf) == true
    @test isnan(x) == false
    @test isnan(x * NaN) == true
    @test isreal(x) == true
    @test isreal(x * (1 + 2im)) == false
    @test signbit(x) == true
    @test signbit(-x) == false
    @test isempty(x) == false
    @test isempty(GenericQuantity([0.0, 1.0])) == false
    @test isempty(GenericQuantity(Float64[])) == true
    @test iseven(Quantity(2, length=1)) == true
    @test iseven(Quantity(3, length=1)) == false
    @test isodd(Quantity(2, length=1)) == false
    @test isodd(Quantity(3, length=1)) == true
    @test isone(Quantity(1, length=1)) == true
    @test isone(Quantity(2, length=1)) == false
    @test isinteger(Quantity(2, length=1)) == true
    @test isinteger(Quantity(2.1, length=1)) == false
    @test ispow2(Quantity(2, length=1)) == true
    @test ispow2(Quantity(3, length=1)) == false

    @test nextfloat(x) == Quantity(nextfloat(-1.2), length=2 // 5)
    @test prevfloat(x) == Quantity(prevfloat(-1.2), length=2 // 5)

    y = Quantity(-1, mass=1)
    @test unsigned(y) == Quantity(unsigned(-1), mass=1)

    z = Quantity(-0.3)
    @test float(z) == z
    @test convert(Float32, z) === convert(Float32, -0.3)

end

@testset "Complex numbers" begin
    x = (0.5 + 0.6im) * u"km/s"
    @test string(x) == "(500.0 + 600.0im) m s⁻¹"
    @test typeof(x) == Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test typeof(x^2) == Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test x^2/u"km/s"^2 ≈ (0.5 + 0.6im)^2
    @test x^2.5 ≈ (-9.896195997465055e6 + 1.38810912834778e7im) * u"m^(5/2) * s^(-5/2)"
    @test isreal(x) == false
    @test isreal(abs2(x)) == true
    @test real(x) == 0.5 * u"km/s"
    @test imag(x) == 0.6 * u"km/s"
    @test conj(x) == (0.5 - 0.6im) * u"km/s"
    @test angle(x) == angle(ustrip(x))
    @test adjoint(ustrip(x^2)) ≈ adjoint(x^2) / u"m/s"^2

    # Can create by division as well:
    x = RealQuantity(1.0u"km/s") / (1.0 + 0.5im)
    @test typeof(x) == Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test ustrip(x) ≈ 1000.0 / (1.0 + 0.5im)
    @test ulength(x) == 1.0
    @test utime(x) == -1.0

    x = (1.0 + 0.5im) / RealQuantity(1.0u"km/s")
    @test typeof(x) == Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test ustrip(x) ≈ (1.0 + 0.5im) / 1000.0
    @test ulength(x) == -1.0
    @test utime(x) == 1.0

    # Can also construct using `complex`
    x = complex(1.0u"km/s", 0.5u"km/s")
    @test typeof(x) === Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test x == (1.0 + 0.5im) * u"km/s"

    x2 = complex(RealQuantity(1.0u"km/s"), RealQuantity(0.5u"km/s"))
    @test typeof(x2) === Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test x2 == (1.0 + 0.5im) * u"km/s"

    x3 = complex(RealQuantity(1.0u"km/s"), GenericQuantity(0.5u"km/s"))
    @test typeof(x3) === GenericQuantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test x3 == (1.0 + 0.5im) * u"km/s"

    # Or, construct from quantity and number
    x4 = complex(Quantity(1.0), 0.5)
    @test typeof(x4) === Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test x4 == (1.0 + 0.5im)

    x5 = complex(1.0, Quantity(0.5))
    @test typeof(x5) === Quantity{Complex{Float64}, DEFAULT_DIM_TYPE}
    @test x5 == (1.0 + 0.5im)

    # Error checking
    @test_throws DimensionError complex(1.0u"km/s", 0.5u"m")
    @test_throws DimensionError complex(1.0u"m", 0.5)
    @test_throws DimensionError complex(1.0, 0.5u"m")
end

@testset "Fallbacks" begin
    @test ustrip(0.5) == 0.5
    @test ustrip(ones(32)) == ones(32)
    @test dimension(Dimensions()) === Dimensions()

    @test_throws MethodError Dimensions(1.0)
    @test_throws ErrorException ustrip(Dimensions())
end

@testset "Arrays" begin
    for T in [Float16, Float32, Float64], R in [Rational{Int16}, Rational{Int32}, SimpleRatio{Int}, SimpleRatio{SafeInt16}]
        D = Dimensions{R}

        X = randn(T, 10)
        uX = X .* Quantity{T,D}(1, length=2.5, luminosity=0.5)

        @test eltype(uX) <: Quantity{T,D}
        @test typeof(sum(uX)) <: Quantity{T,D}
        @test sum(X) == ustrip(sum(uX))
        @test dimension(prod(uX)) == dimension(prod([Quantity(T(1), D, length=2.5, luminosity=0.5) for i in 1:10]))
        @test dimension(prod(uX)) == dimension(prod([Quantity(T(1), D, length=2.5, luminosity=0.5) for i in 1:10]))
        @test typeof(dimension(prod(uX))) <: D

        uX = X .* Quantity(2, length=2.5, luminosity=0.5)
        @test sum(X) == 0.5 * ustrip(sum(uX))

        @test isapprox(uX, uX, atol=Quantity{T,D}(1e-6, length=2.5, luminosity=0.5))
        @test_throws DimensionError isapprox(uX, uX, atol=1e-6)

        x = GenericQuantity(ones(T, 32))
        @test ustrip(x + ones(T, 32))[32] == 2
        @test typeof(x + ones(T, 32)) <: GenericQuantity{Vector{T}}
        @test typeof(x - ones(T, 32)) <: GenericQuantity{Vector{T}}
        @test typeof(ones(T, 32) * GenericQuantity(T(1), D, length=1)) <: GenericQuantity{Vector{T}}
        @test typeof(ones(T, 32) / GenericQuantity(T(1), D, length=1)) <: GenericQuantity{Vector{T}}
        @test ones(T, 32) / GenericQuantity(T(1), length=1) == GenericQuantity(ones(T, 32), length=-1)
    end

    @testset "isapprox" begin
        A = QuantityArray([1, 2, 3], u"m")
        B = QuantityArray([1, 2, 3], u"m")
        @test isapprox(A, B)

        A = QuantityArray([1, 2, 3], u"m")
        B = QuantityArray([1, 2, 3], u"s")
        @test_throws DimensionError isapprox(A, B)

        A = QuantityArray([1, 2, 3], u"m")
        B = QuantityArray([1.001, 2.001, 3.001], u"m")
        @test !isapprox(A, B, atol=1e-4u"m")
        @test isapprox(A, B, atol=1e-2u"m")

        A = QuantityArray([1, 2, 3], u"m")
        B = QuantityArray([1.1, 2.1, 3.1], u"m")
        @test !isapprox(A, B, rtol=0.01)
        @test isapprox(A, B, rtol=0.05)

        A = [1u"m", 2u"m", 3u"m"]
        B = QuantityArray([1, 2, 3], u"m")
        @test isapprox(A, B)

        A = [1u"m", 2u"m", 3u"s"]
        B = QuantityArray([1, 2, 3], u"m")
        @test_throws DimensionError isapprox(A, B)
        @test_throws DimensionError isapprox(B, A)

        A = [1u"m", 2u"m"]
        B = [1u"m", 2u"s"]
        @test_throws DimensionError isapprox(A, B)
        @test_throws DimensionError isapprox(B, A)

        # With different rtoldefault:
        A = QuantityArray([1, 2, 3], Quantity{Float16}(u"m"))
        B = QuantityArray([1.01, 2.01, 3.01], Quantity{Float64}(u"m"))
        @test isapprox(A, B)
        @test isapprox(B, A)  # Because we get it from Float16
        @test !isapprox(Quantity{Float64}.(A), B)
        @test !isapprox(B, Quantity{Float64}.(A))

        # With explicit atol=0, rtol=0
        A = [1u"m", 2u"m", 3u"m"]
        B = [1u"m", 2u"m", 3u"m"]
        @test isapprox(A, B, atol=0u"m", rtol=0)
        B = [1.00000001u"m", 2u"m", 3u"m"]
        @test !isapprox(A, B, atol=0u"m", rtol=0)
    end

    x = randn(32) .* u"km/s"
    @test ustrip.(x) == [ustrip(xi) for xi in x]
    @test dimension.(x) == [dimension(u"km/s") for xi in x]
    @test_throws DimensionError dimension([u"km/s", u"km"])

    @test norm(x, 2) ≈ norm(ustrip.(x), 2) * u"m/s"
    @test norm(GenericQuantity(ustrip.(x), length=1, time=-1), 2) ≈ norm(ustrip.(x), 2) * u"m/s"

    @test ustrip(x') == ustrip(x)'
end

@testset "Ranges" begin
    @testset "Ranges from units" begin
        x = [xi for xi in 0.0u"km/s":0.1u"km/s":1.0u"km/s"]
        @test x[2] == 0.1u"km/s"
        @test x[end] == 1.0u"km/s"

        # https://github.com/JuliaLang/julia/issues/56610
        c = collect(1u"inch":0.25u"inch":4u"inch")
        @test c[1] == 1u"inch"
        @test c[end] <= 4u"inch"

        # Test dimensionless quantities
        x = collect(1u"1":2u"1":5u"1")
        @test x == [1, 3, 5] .* u"1"
        @test eltype(x) <: Quantity

        # Test error for missing step
        @test_throws "must specify a step" 1u"km":2u"km"
        @test_throws "must specify a step" 1us"km":2us"km"

        # However, for backwards compatibility, dimensionless ranges are allowed:
        x = collect(1u"1":5u"1")
        @test x == [1, 2, 3, 4, 5]
        @test eltype(x) <: Quantity{Float64}

        # Test errors for incompatible units
        @test_throws DimensionError 1u"km":1u"s":5u"km"
        @test_throws DimensionError 1u"km":1u"m":5u"s"
        @test_throws DimensionError 1u"km":1u"km/s":5u"km"

        # Same for symbolic units!
        @test_throws DimensionError 1us"km":1us"m":5us"inch"
        @test length(1u"inch":1u"m":5u"km") == 5000

        # Test with symbolic units
        x = collect(1us"inch":0.25us"inch":4us"inch")
        @test x[1] == 1us"inch"
        @test x[2] == 1.25us"inch"
        @test x[end] == 4us"inch"
    end

    @testset "Multiplying ranges with units" begin
        # Test multiplying ranges with units
        x = (1:0.25:4)u"inch"
        @test x isa StepRangeLen
        @test first(x) == 1u"inch"
        @test x[2] == 1.25u"inch"
        @test last(x) == 4u"inch"
        @test length(x) == 13

        # Integer range (but real-valued unit)
        x = (1:4)u"inch"
        @test x isa StepRangeLen
        @test first(x) == 1u"inch"
        @test x[2] == 2u"inch"
        @test last(x) == 4u"inch"
        @test length(x) == 4
        @test collect(x)[3] == 3u"inch"

        # Test with floating point range
        x = (1.0:0.5:3.0)u"m"
        @test x isa StepRangeLen
        @test first(x) == 1.0u"m"
        @test x[2] == 1.5u"m"
        @test last(x) == 3.0u"m"
        @test length(x) == 5
        @test collect(x)[3] == 2.0u"m"

        x = (0:0.1:1)u"m"
        @test length(x) == 11
        @test collect(x)[3] == 0.2u"m"

        # Test with symbolic units
        x = (1:0.25:4)us"inch"
        @test x isa StepRangeLen{<:Quantity{Float64,<:SymbolicDimensions}}
        @test first(x) == us"inch"
        @test x[2] == 1.25us"inch"
        @test last(x) == 4us"inch"
        @test length(x) == 13

        # Test that symbolic units preserve their symbolic nature
        x = (0:0.1:1)us"km/h"
        @test x isa AbstractRange
        @test first(x) == 0us"km/h"
        @test x[2] == 0.1us"km/h"
        @test last(x) == 1us"km/h"
        @test length(x) == 11

        # Similarly, integers should stay integers:
        x = (1:4)us"inch"
        @test_skip x isa StepRangeLen{<:Quantity{Int64,<:SymbolicDimensions}}
        @test first(x) == us"inch"
        @test x[2] == 2us"inch"
        @test last(x) == 4us"inch"
        @test length(x) == 4

        # With RealQuantity:
        @test_skip (1.0:4.0) * RealQuantity(u"inch") isa StepRangeLen{<:RealQuantity{Float64,<:SymbolicDimensions}}
        # TODO: This is not available as TwicePrecision interacts with Real in a way
        #       that demands many other functions to be defined.
    end
end

@testset "Alternate dimension construction" begin
    z = Quantity(-52, length=1) * Dimensions(mass=2)
    z2 = Dimensions(mass=2) * Quantity(-52, length=1)

    @test typeof(z).parameters[1] <: Int
    @test z == z2
    @test ustrip(z) == -52
    @test dimension(z) == Dimensions(length=1, mass=2)
    @test float(z / (z * -1 / 52)) ≈ ustrip(z)

    # Invalid ways to make a quantity:
    @test_throws ErrorException Dimensions(length=1) / 0.5 == Quantity(2.0, length=1)
    @test_throws ErrorException 0.5 / Dimensions(length=1)
    @test_throws ErrorException Dimensions(length=1) * 0.5

    @test 0.5 / Quantity(1, length=1) == Quantity(0.5, length=-1)
    @test 0.5 * Quantity(1, length=1) == Quantity(0.5, length=1)
    @test Quantity(0.5) / Dimensions(length=1) == Quantity(0.5, length=-1)
    @test Quantity(0.5, length=2) / Dimensions(length=1) == Quantity(0.5, length=1)
    @test Dimensions(length=1) / Quantity(0.5, length=2, mass=-5) == Quantity(2, length=-1, mass=5)

    @test Dimensions{Int8}(zeros(Int, 7)...) == Dimensions{Int8}()

    @test zero(Quantity(0.0+0.0im)) + Quantity(1) == Quantity(1.0+0.0im, length=Int8(0))
    @test oneunit(Quantity(0.0+0.0im)) - Quantity(1) == Quantity(0.0+0.0im, length=Int8(0))
    @test typeof(one(Dimensions{Int16})) == Dimensions{Int16}
    @test one(Dimensions{Int16}) == Dimensions(mass=Int16(0))

    @test zero(Quantity(0.0im)) == Quantity(0.0+0.0im)
    @test one(Quantity{ComplexF64}) == Quantity(1.0+0.0im)

    @test zero(Quantity(0.0)) == Quantity(0.0)
    @test typeof(zero(Quantity(0.0))) == Quantity{Float64,DEFAULT_DIM_TYPE}
    @test oneunit(Quantity(1.0)) - Quantity(1.0) == Quantity(0.0)
    @test typeof(one(Quantity(1.0))) == Quantity{DEFAULT_VALUE_TYPE,DEFAULT_DIM_TYPE}
    @test one(Dimensions) == Dimensions()
    @test one(Dimensions()) == Dimensions()
    @test typeof(one(Quantity)) == Quantity{DEFAULT_VALUE_TYPE,DEFAULT_DIM_TYPE}
    @test ustrip(one(Quantity)) === one(DEFAULT_VALUE_TYPE)
    @test typeof(one(GenericQuantity(ones(32, 32)))) == GenericQuantity{Matrix{Float64},DEFAULT_DIM_TYPE}
    @test dimension(one(GenericQuantity(ones(32, 32), length=1))) == Dimensions()

    x = Quantity(1, length=1)

    @test zero(x) == Quantity(0, length=1)
    @test typeof(zero(x)) == Quantity{Int64,DEFAULT_DIM_TYPE}

    # Invalid calls:
    @test_throws ErrorException zero(Quantity)
    @test_throws ErrorException zero(Dimensions())
    @test_throws ErrorException zero(Dimensions)
    @test_throws ErrorException oneunit(Quantity)
    @test_throws ErrorException oneunit(Dimensions())
    @test_throws ErrorException oneunit(Dimensions)

    @test sqrt(z * -1) == Quantity(sqrt(52), length=1 // 2, mass=1)
    @test cbrt(z) == Quantity(cbrt(-52), length=1 // 3, mass=2 // 3)

    @test_throws ErrorException 1.0 * (Dimensions(length=3)^2)

    x = 0.9u"km/s"
    y = 0.3 * x
    @test x > y
    @test y < x

    x = Quantity(1.0)

    @test x == 1.0
    @test x >= 1.0
    @test x < 2.0

    @test_throws DimensionError x < 1.0u"m"
end

@testset "Manual construction" begin
    d = Dimensions(length=-0.2, luminosity=2)
    q = Quantity(0.5, inv(d))
    @test q == Quantity(0.5, length=0.2, luminosity=-2)
end

@testset "Conversions" begin
    d = Dimensions(Rational{Int16}, mass=2)
    d32 = convert(Dimensions{Rational{Int32}}, d)
    @test typeof(d) == Dimensions{Rational{Int16}}
    @test typeof(d32) == Dimensions{Rational{Int32}}
    @test umass(d) == 2
    @test umass(d32) == 2
    @test typeof(umass(d32)) == Rational{Int32}

    # Should not change:
    @test convert(Dimensions, d) === d

    q = Quantity(0.5, d)
    q32_32 = convert(Quantity{Float32,Dimensions{Rational{Int32}}}, q)
    @test typeof(q) == Quantity{Float64,Dimensions{Rational{Int16}}}
    @test typeof(q32_32) == Quantity{Float32,Dimensions{Rational{Int32}}}
    @test ustrip(q) == 0.5
    @test ustrip(q32_32) == 0.5
    @test typeof(ustrip(q)) == Float64
    @test typeof(ustrip(q32_32)) == Float32
    @test dimension(q32_32) == dimension(q)
    @test umass(q) == 2
    @test umass(q32_32) == 2
    @test typeof(umass(q32_32)) == Rational{Int32}
    @test typeof(convert(Quantity{Float16}, q)) == Quantity{Float16,Dimensions{Rational{Int16}}}
    @test convert(Quantity, q) === q

    # Test that regular type promotion applies:
    q = Quantity(2, d)
    @test typeof(q) == Quantity{Int64,typeof(d)}
    @test typeof(q ^ 2) == Quantity{Int64,typeof(d)}
    @test typeof(0.5 * q) == Quantity{Float64,typeof(d)}
    @test typeof(inv(q)) == Quantity{Float64,typeof(d)}

    # Automatic conversions via constructor:
    for T in [Float16, Float32, Float64, BigFloat], R in [DEFAULT_DIM_BASE_TYPE, Rational{Int16}, Rational{Int32}, SimpleRatio{Int}, SimpleRatio{SafeInt16}]
        D = Dimensions{R}
        q = Quantity{T,D}(2, length=1.5)
        @test typeof(q) == Quantity{T,D}
        @test typeof(ustrip(q)) == T
        @test typeof(ulength(q)) == R

        # Now, without R, the default will be DEFAULT_DIM_BASE_TYPE:
        q = Quantity{T}(2, length=1.5)
        @test typeof(q) == Quantity{T,DEFAULT_DIM_TYPE}
        @test typeof(ustrip(q)) == T
        @test typeof(ulength(q)) == DEFAULT_DIM_BASE_TYPE

        # Just dimensions:
        d = D(length=1.5)
        @test typeof(d) == D
        @test typeof(ulength(d)) == R
    end

    # Can convert type with float(...):
    x = Quantity{Int}(u"m/s")
    @test typeof(x) == Quantity{Int64,DEFAULT_DIM_TYPE}
    y = float(typeof(x))
    @test y == Quantity{Float64,DEFAULT_DIM_TYPE}
end

@testset "Units" begin
    x = 1.3u"km/s^2"
    @test ustrip(x) == 1300  # SI base units
    @test ulength(x) == 1
    @test utime(x) == -2

    y = 0.9u"sqrt(mΩ)"
    @test typeof(y) == with_type_parameters(DEFAULT_QUANTITY_TYPE, Float64, DEFAULT_DIM_TYPE)
    @test ustrip(y) ≈ 0.02846049894151541
    @test ucurrent(y) == -1
    @test ulength(y) == 1

    y = BigFloat(0.3) * u"mΩ"
    @test typeof(y) == with_type_parameters(DEFAULT_QUANTITY_TYPE, BigFloat, DEFAULT_DIM_TYPE)
    @test ustrip(y) ≈ 0.0003
    @test ulength(y) == 2

    y32 = convert(with_type_parameters(DEFAULT_QUANTITY_TYPE, Float32, Dimensions{Rational{Int16}}), y)
    @test typeof(y32) == with_type_parameters(DEFAULT_QUANTITY_TYPE, Float32, Dimensions{Rational{Int16}})
    @test ustrip(y32) ≈ 0.0003

    z = u"yr"
    @test utime(z) == 1
    @test ustrip(z) ≈ 60 * 60 * 24 * 365.25
    @test z == uparse("yr")

    # Test type stability of extreme range of units
    @test typeof(u"1") == DEFAULT_QUANTITY_TYPE
    @test typeof(u"1f0") == DEFAULT_QUANTITY_TYPE
    @test typeof(u"s"^2) == DEFAULT_QUANTITY_TYPE
    @test typeof(u"Ω") == DEFAULT_QUANTITY_TYPE
    @test typeof(u"Gyr") == DEFAULT_QUANTITY_TYPE
    @test typeof(u"fm") == DEFAULT_QUANTITY_TYPE
    @test typeof(u"fm"^2) == DEFAULT_QUANTITY_TYPE

    @test_throws ErrorException eval(:(u":x"))

    @test_throws "Symbol x not found" uparse("x")
    @test_throws "Symbol c found in `Constants` but not `Units`" uparse("c")
    @test_throws "Unexpected expression" uparse("import ..Units")
    @test_throws "Unexpected expression" uparse("(m, m)")
    @test_throws LoadError eval(:(us"x"))
    @test_throws "Symbol x not found" sym_uparse("x")
    @test_throws "Symbol c found in `Constants` but not `Units`" sym_uparse("c")
    @test_throws "Unexpected expression" sym_uparse("import ..Units")
    @test_throws "Unexpected expression" sym_uparse("(m, m)")
end

@testset "Constants" begin
    @test Constants.h * Constants.c / (1000.0u"nm") ≈ 1.9864458571489284e-19u"J"

    # Compute period of Earth based on solar mass and semi-major axis:
    a = u"Constants.au"
    @test isapprox(sqrt(4π^2 * a^3 / (Constants.G * Constants.M_sun)), 1u"yr"; rtol=1e-3)
end

@testset "Additional tests of FixedRational" begin
    @test convert(Int64, FixedRational{Int64,1000}(2 // 1)) == 2
    @test convert(Int32, FixedRational{Int64,1000}(3 // 1)) == 3
    @test convert(Bool, FixedRational{Int8,6}(1//1)) === true
    @test convert(Bool, FixedRational{Int8,6}(0//1)) === false

    @test_throws InexactError convert(Int32, FixedRational{Int8,6}(2//3))
    @test_throws InexactError convert(Bool, FixedRational{Int8,6}(2//1))

    @test_throws "Refusing to" promote(FixedRational{Int,10}(2), FixedRational{Int,4}(2))

    f64 = FixedRational{Int,10}(2)
    f8 = FixedRational{Int8,10}(2)
    @test promote(f64, f8) == (2, 2)
    @test typeof(promote(f64, f8)) == typeof((f64, f64))
    @test typeof(promote(FixedRational{Int8,10}(2), FixedRational{Int8,10}(2))) == typeof((f8, f8))
    @test promote_type(Float64, typeof(f64)) == Float64

    # Required to hit integer branch (otherwise will go to `literal_pow`)
    f(i::Int) = Dimensions(length=1, mass=-1)^i
    @test f(2) == Dimensions(length=2, mass=-2)

    # Null conversion
    @test typeof(FixedRational{Int,10}(FixedRational{Int,10}(2))) == FixedRational{Int,10}

    # Conversion to Rational without specifying type
    @test convert(Rational, FixedRational{UInt8,6}(2)) === Rational{UInt8}(2)

    @test record_show(FixedRational{Int,10}(2)) == "2"
    @test record_show(FixedRational{Int,10}(11//10)) == "11//10"

    # Promotion rules
    @test promote_type(FixedRational{Int64,10},FixedRational{BigInt,10}) == FixedRational{BigInt,10}
    @test promote_type(Rational{Int8}, FixedRational{Int,12345}) == Rational{Int}
    @test promote_type(Int8, FixedRational{Int,12345}) == FixedRational{Int,12345}

    # Bug where user would create a FixedRational{::Type{Int32}, ::Int64} and get stack overflow,
    # because the stored type was FixedRational{::Type{Int32}, ::Int32}
    x = 10u"m"
    user_quantity = Quantity(10.0, Dimensions{FixedRational{Int32,25200}}(1, 0, 0, 0, 0, 0, 0))
    @test x == user_quantity
end

@testset "Quantity promotion" begin
    q1 = Quantity(1.0, length=1)
    q2 = Quantity(2, mass=1)
    @test typeof(promote(q1, q2)) == typeof((q1, q1))

    q = 0.5u"km/s"
    x = [0.5, q]
    @test x isa Vector{typeof(q)}
    @test x[1] == convert(typeof(q), 0.5)

    q = GenericQuantity(0.5u"km/s")
    x = [0.5, q]
    @test x isa Vector{typeof(q)}

    # Promotion with custom numeric type:
    @eval struct MyNumber <: Real
        x::Float64
    end
    a = RealQuantity(0.5u"km/s")
    b = MyNumber(0.5)
    ar = [a, b]
    @test ar isa Vector{Real}
    @test a === ar[1]
    @test b === ar[2]
    @test promote_type(MyNumber, typeof(a)) == Real

    # Explicit conversion so coverage can see it:
    D = DEFAULT_DIM_TYPE
    @test promote_type(Quantity{Float32,D}, Float64) == Quantity{Float64,D}
    @test promote_type(Quantity{Float32,D}, Quantity{Float64,D}) == Quantity{Float64,D}
    @test promote_type(Quantity{Float32,D}, GenericQuantity{Float64,D}) == GenericQuantity{Float64,D}
    @test promote_type(GenericQuantity{Float32,D}, GenericQuantity{Float64,D}) == GenericQuantity{Float64,D}
    @test promote_type(SymbolicDimensions{Rational{Int}}, SymbolicDimensions{DEFAULT_DIM_BASE_TYPE}) == SymbolicDimensions{Rational{Int}}
    @test promote_type(Dimensions{Rational{Int}}, SymbolicDimensions{DEFAULT_DIM_BASE_TYPE}) == Dimensions{Rational{Int}}

    @test promote_quantity_on_quantity(RealQuantity, RealQuantity) == RealQuantity
end

struct MyDimensions{R} <: AbstractDimensions{R}
    length::R
    mass::R
    time::R
end
struct MyQuantity{T,D<:AbstractDimensions} <: AbstractQuantity{T,D}
    value::T
    dimensions::D
end

@testset "Custom dimensions" begin
    for T in [Float32, Float64], R in [Rational{Int64}, Rational{Int32}]
        D = MyDimensions{R}
        x = MyQuantity(T(0.1), D, length=0.5)
        @test x * x ≈ MyQuantity(T(0.01), D, length=1)
        @test typeof(x * x) == MyQuantity{T,D}
        @test one(MyQuantity{T,D}) == MyQuantity(one(T), MyDimensions(R))
        @test zero(x) == MyQuantity(zero(T), D, length=0.5)
        @test oneunit(x) + x == MyQuantity(T(1.1), D, length=0.5)
        @test typeof(oneunit(x) + x) == MyQuantity{T,D}

        # Automatic conversions:
        @test typeof(MyQuantity(1, MyDimensions, length=0.5)) == MyQuantity{typeof(1),MyDimensions{DEFAULT_DIM_BASE_TYPE}}
        @test typeof(MyQuantity{T}(1, MyDimensions, length=0.5)) == MyQuantity{T,MyDimensions{DEFAULT_DIM_BASE_TYPE}}
        @test typeof(MyQuantity{T}(1, D, length=0.5)) == MyQuantity{T,D}
        @test typeof(MyQuantity{T,D}(0.1, length=0.5)) == MyQuantity{T,D}
        @test typeof(0.5 * MyQuantity{T,D}(0.1, length=0.5)) == MyQuantity{promote_type(T,Float64),D}

        x = MyQuantity(big(0.1), length=1)
        @test typeof(x) == MyQuantity{BigFloat,DEFAULT_DIM_TYPE}
        @test typeof(MyQuantity{T}(x)) == MyQuantity{T,DEFAULT_DIM_TYPE}

        # Using MyDimensions inside regular Quantity:
        x = Quantity(T(0.1), MyDimensions(R, length=0.5))
        @test typeof(x) == Quantity{T,MyDimensions{R}}
        @test typeof(x * x) == Quantity{T,MyDimensions{R}}
        @test ulength(x * x) == 1
        @test dimension(x * x) == MyDimensions(R, length=1)

        # Errors:
        @test_throws ErrorException zero(MyQuantity{T,D})
        @test_throws ErrorException oneunit(MyQuantity{T,D})
        @test_throws ErrorException 1.0 * MyDimensions()
    end
    @test MyQuantity(0.1, DEFAULT_DIM_TYPE, length=0.5) == MyQuantity(0.1, length=0.5)
    @test MyQuantity(0.1, DEFAULT_DIM_TYPE, length=0.5) == MyQuantity(0.1, length=1//2)

    # Can construct using args directly:
    @test typeof(MyDimensions(1, 1, 1)) == MyDimensions{Int}
    @test typeof(MyDimensions{Float64}(1, 1, 1)) == MyDimensions{Float64}

    # Can use the default constructorof, with_type_parameters, and dimension_names:
    @test DynamicQuantities.constructorof(MyDimensions{Float64}) == MyDimensions
    @test DynamicQuantities.constructorof(MyQuantity{Float64}) == MyQuantity
    @test DynamicQuantities.with_type_parameters(MyDimensions{Float64}, Rational{Int}) == MyDimensions{Rational{Int}}
    @test DynamicQuantities.with_type_parameters(MyQuantity{Float64,DEFAULT_DIM_TYPE}, Float32, MyDimensions{Float64}) == MyQuantity{Float32,MyDimensions{Float64}}
    @test DynamicQuantities.dimension_names(MyDimensions{Float64}) == (:length, :mass, :time)

    # But, we always need to use a quantity when mixing with mathematical operations:
    @test_throws ErrorException MyQuantity(0.1) + 0.1 * MyDimensions()

    # Explicitly test that `promote_quantity_on_quantity` has a reasonable default
    @test promote_quantity_on_quantity(typeof(MyQuantity(0.1)), typeof(MyQuantity(0.1))) == MyQuantity{Float64,DEFAULT_DIM_TYPE}
end

@testset "Symbolic dimensions" begin
    # TODO: Remove constructors for sym3 and sym4s?
    sym1 = @inferred(SymbolicDimensions(; m=3, s=-1))
    sym2 = @inferred(SymbolicDimensions{Rational{Int}}(; m=3, s=-1))
    sym3 = @inferred(SymbolicDimensions(Rational{Int}; m=3, s=-1))
    sym4 = @inferred(SymbolicDimensions{Int}(Rational{Int}; m=3, s=-1))
    for (sym, T) in (
        (sym1, DynamicQuantities.DEFAULT_DIM_BASE_TYPE), (sym2, Rational{Int}), (sym3, Rational{Int}), (sym4, Rational{Int}),
    )
        @test sym isa SymbolicDimensions{T}

        # Properties
        @test sym.m == 3
        @test sym.s == -1
        @test propertynames(sym) == DynamicQuantities.ALL_SYMBOLS
        @test issubset((:m, :s), propertynames(sym))
        @test all(propertynames(sym)) do x
            val = getproperty(sym, x)
            return x === :m ? val == 3 : (x === :s ? val == -1 : iszero(val))
        end

        # Internal constructor
        @test DynamicQuantities.constructorof(typeof(sym)) === SymbolicDimensions

        # Equality comparisons
        @test sym == sym
        @test sym == copy(sym)
        @test sym !== copy(sym)
        @test sym == SymbolicDimensions{Int}(; m=3, s=-1)
        @test SymbolicDimensions{Int}(; m=3, s=-1) == sym
        @test sym == SymbolicDimensions(; m=3, g=0, s=-1)
        @test SymbolicDimensions(; m=3, g=0, s=-1) == sym
        @test sym == SymbolicDimensions(; m=3, s=-1, K=0)
        @test SymbolicDimensions(; m=3, s=-1, K=0) == sym
        @test sym != SymbolicDimensions(; m=2, s=-1)
        @test SymbolicDimensions(; m=2, s=-1) != sym
        @test sym != SymbolicDimensions(; m=3, g=1, s=-1)
        @test SymbolicDimensions(; m=3, g=1, s=-1) != sym
        @test sym != SymbolicDimensions(; m=3, s=-1, K=1)
        @test SymbolicDimensions(; m=3, s=-1, K=1) != sym

        @test !iszero(sym)
    end

    q = 1.5us"km/s"
    @test q == 1.5 * us"km" / us"s"
    @test typeof(q) <: with_type_parameters(DEFAULT_QUANTITY_TYPE, Float64, SymbolicDimensions{DEFAULT_DIM_BASE_TYPE})
    @test string(dimension(q)) == "s⁻¹ km"
    @test uexpand(q) == 1.5u"km/s"
    @test string(dimension(us"Constants.au^1.5")) == "au³ᐟ²"
    @test string(dimension(uexpand(us"Constants.au^1.5"))) == "m³ᐟ²"
    @test uexpand(2.3us"Constants.au^1.5") ≈ 2.3u"Constants.au^1.5"
    @test iszero(dimension(us"1.0")) == true
    @test uexpand(inv(us"Constants.au")) ≈ 1/u"Constants.au"
    @test dimension(inv(us"s") * us"km") == dimension(us"km/s")
    @test dimension(inv(us"s") * us"m") != dimension(us"km/s")
    @test dimension(uexpand(inv(us"s") * us"m")) == dimension(uexpand(us"km/s"))

    f2(i::Int) = us"s"^i
    @inferred f2(5)
    @test uexpand(f2(5)) == u"s"^5

    @test_throws ErrorException uparse("'c'")
    @test_throws ErrorException sym_uparse("'c'")

    # For constants which have a namespace collision, the numerical expansion is used:
    @test dimension(us"Constants.au")[:au] == 1
    @test dimension(us"Constants.h")[:h] == 0
    @test dimension(us"h")[:h] == 1

    @test us"Constants.h" != us"h"
    @test uexpand(us"Constants.h") == u"Constants.h"

    # Actually expands to:
    @test string(dimension(us"Constants.h")) == "h_constant"

    # Constants have different numerical value from non-symbolic one:
    @test ustrip(us"Constants.h") != ustrip(u"Constants.h")
    @test ustrip(us"Constants.au") != ustrip(u"Constants.au")

    # Test conversion
    @test typeof(SymbolicDimensions{Rational{Int}}(dimension(us"km/s"))) == SymbolicDimensions{Rational{Int}}
    @test convert(Quantity{Float64,SymbolicDimensions}, u"kg") == 1.0us"kg"
    @test convert(Quantity{Float64,SymbolicDimensions}, u"cm") == 1e-2us"m"
    @test convert(Quantity{Float64,Dimensions}, 3.5us"kg/s") == 3.5u"kg/s"
    @test convert(Quantity{Float64,Dimensions}, 3.5us"Constants.pc") == 3.5u"Constants.pc"

    # Helpful error if symbol not found:
    sym5 = dimension(us"km/s")
    @test_throws "my_special_symbol is not available as a symbol" sym5.my_special_symbol

    # Test deprecated method
    q = 1.5us"km/s"
    @test expand_units(q) == uexpand(q)

    # Test promotions:
    x = Quantity{Float32,SymbolicDimensions{Rational{Int}}}(0.2us"km/s")
    y = 0.5us"km/s"
    qa = [x, y]
    @test qa isa Vector{Quantity{Float64,SymbolicDimensions{Rational{Int}}}}
    DynamicQuantities.with_type_parameters(SymbolicDimensions{Float64}, Rational{Int}) == SymbolicDimensions{Rational{Int}}

    # Many symbols in one:
    x = us"pm * fm * nm * μm * mm * cm * dm * m * km"
    y = us"s"
    @test inv(x) != x
    @test dimension(inv(x)).pm == -1
    @test x != y
    @test y != x
    @test dimension(uexpand(x * y)) == dimension(u"m^9 * s")
    z = uexpand(x)
    @test x == z

    # Trigger part of map_dimensions missed elsewhere
    x = us"km"
    y = us"km"
    @test x * y |> uexpand == u"km^2"
    @test x / y |> uexpand == u"1"
    @test map_dimensions(+, dimension(us"km"), dimension(us"km")) == dimension(us"km^2")
    @test map_dimensions(-, dimension(us"km"), dimension(us"km")) == dimension(us"1")

    @testset "Promotion with Dimensions" begin
        x = 0.5u"cm"
        y = -0.03u"m"
        x_s = 0.5us"cm"
        for op in (+, -, *, /, atan, atand, copysign, flipsign, mod)
            @test op(x, y) == op(x_s, y)
            @test op(y, x) == op(y, x_s)
        end
    end
end

@testset "uconvert" begin
    @test uconvert(us"nm", 5e-9u"m") ≈ (5e-9u"m" |> uconvert(us"nm")) ≈ 5us"nm"
    @test_throws DimensionError uconvert(us"nm * J", 5e-9u"m")

    # Types:
    @test typeof(uconvert(us"nm", 5e-9u"m")) <: constructorof(DEFAULT_QUANTITY_TYPE){Float64,<:SymbolicDimensions}
    @test typeof(uconvert(us"nm", GenericQuantity(5e-9u"m"))) <: GenericQuantity{Float64,<:SymbolicDimensions}
    @test uconvert(GenericQuantity(us"nm"), GenericQuantity(5e-9u"m")) ≈ 5us"nm"
    @test uconvert(GenericQuantity(us"nm"), GenericQuantity(5e-9u"m")) ≈ GenericQuantity(5us"nm")

    # We only want to convert the dimensions, and ignore the quantity type:
    @test typeof(uconvert(GenericQuantity(us"nm"), 5e-9u"m")) <: constructorof(DEFAULT_QUANTITY_TYPE){Float64,<:SymbolicDimensions}

    q = 1.5u"Constants.M_sun"
    qs = uconvert(us"Constants.M_sun", 5.0 * q)
    @test qs ≈ 7.5us"Constants.M_sun"
    @test dimension(qs)[:kg] == 0
    @test dimension(qs)[:g] == 0
    @test dimension(qs)[:M_sun] == 1
    @test uexpand(qs) ≈ 5.0 * q

    @test dimension(1u"m" |> us"nm")[:nm] == 1
    @test dimension(1u"m" |> us"nm")[:m] == 0

    # Refuses to convert to non-unit quantities:
    @test_throws AssertionError uconvert(1.2us"m", 1.0u"m")
    @test_throws "You passed a quantity" uconvert(1.2us"m", 1.0u"m")

    # Refuses to convert to `Dimensions`:
    @test_throws ErrorException uconvert(1u"m", 5.0us"m")
    @test_throws "You can only `uconvert`" uconvert(1u"m", 5.0us"m")

    for Q in (RealQuantity, Quantity, GenericQuantity)
        # Different types require converting both arguments:
        q = convert(Q{Float16}, 1.5u"g")
        qs = uconvert(convert(Q{Float16}, us"g"), 5 * q)
        @test typeof(qs) <: Q{Float16,<:SymbolicDimensions{<:Any}}
        @test isapprox(qs, 7.5us"g"; atol=0.01us"g")

        # Arrays
        x = [1.0, 2.0, 3.0] .* Q(u"kg")
        xs = x .|> uconvert(us"g")
        @test typeof(xs) <: Vector{<:Q{Float64,<:SymbolicDimensions{<:Any}}}
        @test xs[2] ≈ Q(2000us"g")

        # Arrays
        x2 = [1.0, 2.0, 3.0] .* Q(u"kg")
        xs2 = x2 .|> us"g"
        @test typeof(xs2) <: Vector{<:Q{Float64,<:SymbolicDimensions{<:Any}}}
        @test xs2[2] ≈ Q(2000us"g")
        @test ustrip(xs2[2]) ≈ 2000
        
        x_qa = QuantityArray(x)
        xs_qa = x_qa .|> uconvert(us"g")
        @test typeof(xs_qa) <: QuantityArray{Float64,1,<:SymbolicDimensions{<:Any}}
        @test xs_qa[2] ≈ Q(2000us"g")
        @test ustrip(xs_qa[1]) ≈ 1000

        x_qa1 = QuantityArray(x)
        xs_qa1 = x_qa1 .|> us"g"
        @test typeof(xs_qa1) <: QuantityArray{Float64,1,<:SymbolicDimensions{<:Any}}
        @test xs_qa1[2] ≈ Q(2000us"g")
        @test ustrip(xs_qa1[3]) ≈ 3000

        # Without vectorized call:
        xs_qa2 = x_qa |> uconvert(us"g")
        @test typeof(xs_qa2) <: QuantityArray{Float64,1,<:SymbolicDimensions{<:Any}}
        @test xs_qa2[2] ≈ Q(2000us"g")
        @test ustrip(xs_qa2[2]) ≈ 2000

        xs_qa3 = x_qa |> us"g"
        @test typeof(xs_qa3) <: QuantityArray{Float64,1,<:SymbolicDimensions{<:Any}}
        @test xs_qa3[2] ≈ Q(2000us"g")
        @test ustrip(xs_qa3[3]) ≈ 3000
    end
end

@testset "Test missing" begin
    x = 1.0u"m"
    @test round(Union{Int32,Missing}, FixedRational{Int32,100}(1)) isa Int32
    @test round(Union{Int32,Missing}, FixedRational{Int32,100}(1)) === Int32(1)
    y = missing
    @test isless(x, y) === missing
    @test isless(y, x) === missing
    @test (x == y) === missing
    @test (y == x) === missing
    @test isapprox(x, y) === missing
    @test isapprox(y, x) === missing

    x = 1.0u"m"
    y = missing
    @test x * y === missing
    @test x + y === missing
end


@testset "Test ambiguities" begin
    @testset "FixedRational" begin
        R = DEFAULT_DIM_BASE_TYPE
        x = convert(R, 10)
        y = convert(R, 5)
        @test promote(x, y) == (x, y)
        @test_throws ErrorException promote(x, convert(FixedRational{Int32,100}, 10))
        @test promote_type(typeof(u"km/s"), typeof(convert(Quantity{Float32}, u"km/s"))) <: Quantity{Float64}

        x = FixedRational{Int32,100}(1)
        # Need explicit `promote_rule` calls here so coverage picks it up
        @test promote_rule(typeof(x), typeof(true)) == typeof(x)
        @test promote_rule(typeof(true), typeof(x)) == typeof(x)
        @test promote_rule(typeof(x), typeof(BigFloat(1))) == promote_type(Rational{Int32}, BigFloat)
        @test promote_rule(typeof(BigFloat(1)), typeof(x)) == promote_type(Rational{Int32}, BigFloat)
        @test promote_rule(typeof(x), typeof(π)) == promote_type(Rational{Int32}, typeof(π))
        @test promote_rule(typeof(π), typeof(x)) == promote_type(Rational{Int32}, typeof(π))
    end

    @testset "Unimplemented on purpose" begin
        x = 1.0u"m"
        s = "test"
        y = WeakRef(s)
        @test_throws ErrorException x == y
        @test_throws ErrorException y == x

        @test_throws ErrorException DynamicQuantities.SymbolicDimensionsSingleton{Int}(Dimensions())
        @test_throws ErrorException DynamicQuantities.SymbolicDimensionsSingleton{Int}(UInt32)
    end

    @testset "Arrays" begin
        qarr1 = QuantityArray(randn(3), u"km/s")
        qarr2 = qarr1
        @test convert(typeof(qarr2), qarr2) === qarr1
    end

    @testset "Rational power law" begin
        x = RealQuantity(1.0u"m")
        y = x ^ (3//2)
        @test y == Quantity(1.0, length=3//2)
        @test typeof(y) == RealQuantity{Float64,DEFAULT_DIM_TYPE}
    end

    @testset "Numeric promotion rules" begin
        for Q in (RealQuantity, Quantity, GenericQuantity)
            x = Q(1.0u"m")
            @test promote_type(typeof(x), Bool) == typeof(x)
            @test promote_type(Bool, typeof(x)) == typeof(x)
            @test promote_type(typeof(x), BigFloat) == with_type_parameters(Q, BigFloat, DEFAULT_DIM_TYPE)
            @test promote_type(BigFloat, typeof(x)) == with_type_parameters(Q, BigFloat, DEFAULT_DIM_TYPE)
        end
    end

    @testset "Complex numbers" begin
        for Q in (RealQuantity, Quantity, GenericQuantity)
            x = 1.0im
            y = Q(0.5u"m")
            @test typeof(x * y) == with_type_parameters(promote_quantity_on_value(Q, ComplexF64), Complex{Float64}, DEFAULT_DIM_TYPE)
            @test typeof(y * x) == with_type_parameters(promote_quantity_on_value(Q, ComplexF64), Complex{Float64}, DEFAULT_DIM_TYPE)
            @test ustrip(x * y) == 0.5im
            @test ustrip(y * x) == 0.5im

            # Bool version
            x = true * im
            y = Q(0.5u"m")
            @test typeof(x * y) == with_type_parameters(promote_quantity_on_value(Q, ComplexF64), Complex{Float64}, DEFAULT_DIM_TYPE)
            @test typeof(y * x) == with_type_parameters(promote_quantity_on_value(Q, ComplexF64), Complex{Float64}, DEFAULT_DIM_TYPE)
            @test ustrip(x * y) == 0.5im
            @test ustrip(y * x) == 0.5im

            # Complex powers
            x = Q(0.5u"1")
            out = x ^ (1 + 2im)
            @test typeof(out) == with_type_parameters(promote_quantity_on_value(Q, ComplexF64), Complex{Float64}, DEFAULT_DIM_TYPE)
            @test ustrip(out) ≈ 0.5 ^ (1 + 2im)

            for CT in (Complex, Complex{Bool})
                x = Q(1.0)
                @test CT(x) == CT(1.0)
                @test typeof(CT(x)) <: CT
                x = Q(1.0, length=1)
                @test_throws AssertionError CT(x)
            end
        end
    end

    @testset "Bool" begin
        for Q in (RealQuantity, Quantity, GenericQuantity)
            x = Q(1.0u"1")
            @test Bool(x) == true
            @test Bool(ustrip(x)) == true
            @test Bool(Q(0.0u"1")) == false
            @test Bool(ustrip(Q(0.0u"1"))) == false
            x = Q(1.0u"m")
            @test_throws AssertionError Bool(x)
        end
    end
end

@testset "Arrays" begin
    for Q in (RealQuantity, Quantity, GenericQuantity)
        @testset "Basics $Q" begin
            x = QuantityArray(randn(32), Q(u"km/s"))
            @test ustrip(sum(x)) ≈ sum(ustrip(x))

            # Setting index with different quantity:
            x[5] = Q(5, length=1, time=-1)
            @test x[5] == Q(5, length=1, time=-1)

            y = randn(32)
            @test ustrip(QuantityArray(y, Q(u"m"))) ≈ y

            f_square(v) = v^2 * 1.5 - v^2
            @test sum(f_square.(QuantityArray(y, Q(u"m")))) ≈ sum(f_square.(y) .* Q(u"m^2"))

            y_q = QuantityArray(y, Q(u"m * cd / s"))
            @test typeof(f_square.(y_q)) == typeof(y_q)

            for get_u in (ulength, umass, utime, ucurrent, utemperature, uluminosity, uamount)
                @test get_u(f_square.(y_q)) == get_u(y_q) * 2
                @test get_u(f_square(first(y_q))) == get_u(y_q) * 2
            end

            # Test default constructors:
            @test QuantityArray(ones(3), u"m/s") == QuantityArray(ones(3), length=1, time=-1)
            @test typeof(QuantityArray(ones(3), u"m/s")) <: QuantityArray{Float64,1,<:Dimensions,<:constructorof(DEFAULT_QUANTITY_TYPE),<:Array}

            # We can create quantity arrays with generic quantity
            @test typeof(QuantityArray([[1.0], [2.0, 3.0]], dimension(u"m/s"))) <: QuantityArray{<:Any,1,<:Dimensions,<:GenericQuantity,<:Array}
            @test typeof(QuantityArray([[1.0], [2.0, 3.0]], GenericQuantity(u"m/s"))) <: QuantityArray{<:Any,1,<:Dimensions,<:GenericQuantity,<:Array}
            @test QuantityArray([[1.0], [2.0, 3.0]], GenericQuantity(u"km/s"))[1][1] == 1000u"m/s"

            fv_square(v) = f_square.(v)
            @inferred fv_square(y_q)

            # Same array type:
            s_x = QuantityArray(SArray{Tuple{32}}(ustrip(x)), dimension(x), Q)
            @test typeof(s_x) <: QuantityArray{Float64,1,<:Dimensions,<:Q{Float64,<:Dimensions},<:SArray}
            output_s_x = (xi -> xi^2).(s_x)
            @test array_type(output_s_x) <: MArray
            @test dim_type(output_s_x) <: Dimensions{<:FixedRational}
            @test value_type(output_s_x) == Float64
            @test dimension(output_s_x) == dimension(x)^2
            fv_square2(x) = (xi -> xi^2).(x)
            @inferred fv_square2(s_x)

            # isapprox for QuantityArray's
            @test isapprox(x, x, atol=Q(1e-6u"km/s"))
            @test_throws DimensionError isapprox(x, x, atol=1e-6)
        end

        @testset "Copying $Q" begin
            x = QuantityArray(randn(3), Q(u"km/s"))
            xc = copy(x)
            @test x == xc
            xc[2] *= 0.5
            @test x != xc
        end

        @testset "Utilities $Q" begin
            @test fill(Q(u"m/s"), 10) == QuantityArray(fill(1.0, 10) .* Q(u"m/s"))
            @test ndims(fill(Q(u"m/s"), ())) == 0
            @test fill(Q(u"m/s"), ())[begin] == Q(u"m/s")
        end

        @testset "similar $Q" begin
            qa = QuantityArray(rand(3, 4), Q(u"m"))

            new_qa = similar(qa)
            @test size(new_qa) == size(qa)
            @test eltype(new_qa) == eltype(qa)
            @test dim_type(new_qa) == dim_type(qa)
            @test quantity_type(new_qa) == eltype(new_qa)
            @test dimension(new_qa) == dimension(qa)
            @test isa(ustrip(new_qa), Array{Float64,2})
            @test !isequal(ustrip(qa), ustrip(new_qa))

            new_qa = similar(qa, Float32)
            @test eltype(new_qa) <: Q{Float32}
            @test dim_type(new_qa) == dim_type(qa)
            @test dimension(new_qa) == dimension(qa)
            @test isa(ustrip(new_qa), Array{Float32,2})

            if Q !== GenericQuantity
                new_qa = similar(qa, typeof(GenericQuantity{Float16}(u"km/s")))
                @test eltype(new_qa) <: GenericQuantity{Float16}
                @test dim_type(new_qa) == dim_type(qa)
                @test dimension(new_qa) == dimension(qa)
                @test isa(ustrip(new_qa), Array{Float16,2})
            end

            new_qa = similar(qa, axes(ones(6, 8)))
            @test size(new_qa) == (6, 8)
            @test eltype(new_qa) <: Q{Float64}
            @test dim_type(new_qa) == dim_type(qa)
            @test dimension(new_qa) == dimension(qa)
            @test isa(ustrip(new_qa), Array{Float64,2})

            new_qa = similar(qa, Float32, axes(ones(6, 8)))
            @test size(new_qa) == (6, 8)
            @test eltype(new_qa) <: Q{Float32}

            new_qa = similar(qa, Float32, (6,))
            @test size(new_qa) == (6,)
            @test eltype(new_qa) <: Q{Float32}

            new_qa = similar(qa, (6,))
            @test size(new_qa) == (6,)
            @test eltype(new_qa) <: Q{Float64}

            new_qa = similar(qa, Float32, (6, UInt(3)))
            @test size(new_qa) == (6, 3)
            @test eltype(new_qa) <: Q{Float32}

            new_qa = similar(qa, (6, UInt(3)))
            @test size(new_qa) == (6, 3)
            @test eltype(new_qa) <: Q{Float64}
        end

        @testset "Zero $Q" begin
            qa = QuantityArray(rand(3, 4), Q(u"m"))
            @test zero(qa) == QuantityArray(zeros(3, 4), Q(u"m"))
        end

        @testset "Promotion $Q" begin
            qarr1 = QuantityArray(randn(32), convert(Dimensions{Rational{Int32}}, dimension(u"km/s")), Q)
            qarr2 = QuantityArray(randn(Float16, 32), convert(Dimensions{Rational{Int64}}, dimension(u"km/s")), Q)

            expected_T = Float64
            expected_D = Dimensions{Rational{Int64}}
            expected_type = QuantityArray{expected_T,1,expected_D,Q{Float64,expected_D},Array{expected_T,1}}

            @test promote_type(typeof(qarr1), typeof(qarr2)) == expected_type
            @test typeof(promote(qarr1, qarr2)) == Tuple{expected_type, expected_type}
        end

        @testset "Array concatenation $Q" begin
            qarr1 = QuantityArray(randn(3) .* Q(u"km/s"))
            qarr2 = QuantityArray(randn(3) .* Q(u"km/s"))

            @test ustrip.(hcat(qarr1, qarr2)) == hcat(ustrip(qarr1), ustrip(qarr2))
            @test ustrip.(vcat(qarr1, qarr2)) == vcat(ustrip(qarr1), ustrip(qarr2))
            @test ustrip.(cat(qarr1, qarr2, dims=2)) == cat(ustrip(qarr1), ustrip(qarr2), dims=Val(2))
            @test dimension(hcat(qarr1, qarr2)) == dimension(u"km/s")

            # type stability:
            @inferred hcat(qarr1, qarr2)
            @inferred vcat(qarr1, qarr2)
            @inferred cat(qarr1, qarr2, dims=Val(2))

            # same array type:
            s_qarr1 = QuantityArray(SArray{Tuple{3}}(ustrip(qarr1)), dimension(qarr1), Q)
            s_qarr2 = QuantityArray(SArray{Tuple{3}}(ustrip(qarr2)), dimension(qarr2), Q)
            @test array_type(hcat(s_qarr1, s_qarr2)) <: SArray

            # Test concatenating different arrays:

        end

        @testset "Generic literal_pow $Q" begin
            y = randn(32)
            y_q = QuantityArray(y, Q(u"m"))

            f4(v) = v^4 * 0.3
            @test sum(f4.(QuantityArray(y, Q(u"m")))) ≈ sum(f4.(y) .* Q(u"m^4"))

            f4v(v) = f4.(v)
            @inferred f4v(y_q)
        end

        @testset "Broadcast with single number $Q" begin
            ar1 = QuantityArray(randn(3), Q(u"km/s"))
            @test ustrip(ar1 .* Q(u"m/s")) == ustrip(ar1)
            @test dimension(ar1 .* Q(u"m/s")) == dimension(u"m^2/s^2")
        end

        @testset "Multiple arrays $Q" begin
            ar1 = QuantityArray(randn(3), Q(u"km/s"))
            ar2 = QuantityArray(randn(3, 1), Q(u"km/s"))
            ar3 = randn(3)
            f(x, y, z) = x + y * z
            g(x, y, z) = f.(x, y, z)
            @inferred g(ar1, ar2[:, 1], ar3)
            @test g(ar1, ar2[:, 1], ar3) == [f(ar1[i], ar2[i, 1], ar3[i]) for i in eachindex(ar1)]

            tuple_array = Q.((1u"km/s", 2u"m/s", 3u"cm/s"))
            @inferred g(ar1, tuple_array, ar3)
            @test g(ar1, tuple_array, ar3) == [f(ar1[i], tuple_array[i], ar3[i]) for i in eachindex(ar1)]

            array_of_quantities = Q.([1u"km/s", 2u"m/s", 3u"cm/s"])
            @inferred g(ar1, array_of_quantities, ar3)
            @test g(ar1, tuple_array, ar3) == g(ar1, array_of_quantities, ar3)

            @test typeof(Base.broadcasted(f, ar1, array_of_quantities, 1.0).args[end]) == Float64
            q = Q(u"1")
            @test typeof(Base.broadcasted(f, ar1, array_of_quantities, q).args[end]) == typeof(q)

            # TODO: Type inference here needs to be fixed
            @test_skip @inferred g(ar1, array_of_quantities, Q(u"1"))
            @test g(ar1, array_of_quantities, Q(u"1")) == [f(ar1[i], array_of_quantities[i], 1) for i in eachindex(ar1)]
        end

        @testset "Broadcast nd-arrays $Q" begin
            x = QuantityArray(randn(3, 3), Q(u"A"))
            y = QuantityArray(randn(3, 3), Q(u"cd"))
            @test ustrip(x .* y) == ustrip(x) .* ustrip(y)
        end

        @testset "map $Q" begin
            qa = fill(Q(2.0u"m"), 3)
            @test map(x -> x^2, qa) == QuantityArray(fill(4.0, 3), Q(u"m^2"))

            @test mapreduce(x -> x^2, +, qa) == 12.0u"m^2"
            @inferred mapreduce(x -> x^2, +, qa)
            @test prod(qa) == 8.0u"m^3"
            @inferred prod(qa)

            # Map to non-quantity output:
            @test map(x -> ustrip(x), qa) == fill(2.0, 3)
            @test map(x -> cos(x/dimension(x)), qa) == fill(cos(2.0), 3)

            # Test that we can use a function that returns a different type
            if Q === RealQuantity
                qa = fill(RealQuantity(2.0u"m"), 3)
                @test typeof(qa) <: QuantityArray{Float64,1,<:Dimensions,<:RealQuantity{Float64,<:Dimensions},<:Array}
                qa_complex = map(x -> 1im * x, qa)
                @inferred map(x -> 1im * x, qa)
                @test typeof(qa_complex) <: QuantityArray{ComplexF64,1,<:Dimensions,<:Quantity{ComplexF64,<:Dimensions},<:Array}
                @test dimension(qa_complex) == dimension(u"m")
                @test sum(qa_complex) == 6im * u"m"
            end
        end

        @testset "Collection operations" begin
            # push!, pushfirst!, insert!, append!, prepend!, pop!,
            # popfirst!, deleteat!, popat!, resize!, empty!, sizehint!

            # Test for push!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            push!(arr, 4u"m")
            @test length(arr) == 4
            @test arr[end] == 4u"m"
            @test_throws DimensionError push!(arr, 5u"s")

            # Multiple value
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            push!(arr, 4u"m", 5u"m")
            @test arr == QuantityArray([1u"m", 2u"m", 3u"m", 4u"m", 5u"m"])

            @test_throws DimensionError push!(arr, 4u"m", 5u"m/s")

            # Test for pushfirst!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            pushfirst!(arr, 0u"m")
            @test length(arr) == 4
            @test arr[1] == 0u"m"
            @test_throws DimensionError pushfirst!(arr, 1u"s")

            # Symbolic dimensions also work
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            pushfirst!(arr, 1us"cm")
            @test arr[1] == 1u"cm"

            # Test for insert!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            insert!(arr, 2, 1.5u"m")
            @test length(arr) == 4
            @test arr[2] == 1.5u"m"
            @test_throws DimensionError insert!(arr, 2, 2u"s")

            # Symbolic dimensions
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            insert!(arr, 2, 1.5us"cm")
            @test arr[2] == 1.5u"cm"

            # But not an incompatible symbolic dimensions:
            @test_throws DimensionError insert!(arr, 2, 1.5us"s")

            # Test for append!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            append!(arr, QuantityArray([5u"m", 6u"m"]))
            @test length(arr) == 5
            @test arr[end] == 6u"m"
            @test_throws DimensionError append!(arr, QuantityArray([7u"s", 8u"s"]))

            # Test appending regular array of quantities:
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            append!(arr, [5u"m", 6u"m"])
            @test length(arr) == 5
            @test arr[end] == 6u"m"
            @test_throws DimensionError append!(arr, [7u"s", 8u"s"])

            # Can also append a QuantityArray of symbolic units:
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            orig_type = typeof(arr)
            append!(arr, QuantityArray([5us"cm", 6us"cm"]))
            @test ustrip(arr[end]) == ustrip(6u"cm")
            @test typeof(arr) === orig_type

            # Or, regular units to symbolic units:
            arr = QuantityArray([1us"m", 2us"m", 3us"m"])
            append!(arr, [5u"m", 6u"m"])
            @test dimension(5us"m") == dimension(arr[end-1])

            # Test for prepend!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            prepend!(arr, QuantityArray([-1u"m", -2u"m"]))
            @test length(arr) == 5
            @test arr[1] == -1u"m"
            @test_throws DimensionError prepend!(arr, QuantityArray([-3u"s", -4u"s"]))

            # Test for pop!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            val = pop!(arr)
            @test val == 3u"m"
            @test length(arr) == 2

            # Test for popfirst!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            val = popfirst!(arr)
            @test val == 1u"m"
            @test length(arr) == 2

            # Test for deleteat!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            deleteat!(arr, 2)
            @test length(arr) == 2
            @test arr[2] == 3u"m"

            # Test for popat!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            val = popat!(arr, 2)
            @test val == 2u"m"
            @test length(arr) == 2
            @test arr[1] == 1u"m"
            @test arr[2] == 3u"m"

            # Test for resize!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            resize!(arr, 5)
            @test length(arr) == 5
            @test dimension(arr[end]) == dimension(u"m")

            # Test for empty!
            arr = QuantityArray([1u"m", 2u"m", 3u"m"])
            empty!(arr)
            @test isempty(arr)
            # Still stores the dimension:
            @test dimension(arr) == dimension(u"m")

            # Test for sizehint!
            # There is no easy way to test whether it actually ran,
            # so we create a fake array type that has a custom `sizehint!`
            # which tells us it actually ran.
            isdefined(@__MODULE__, :MyCustomArray) || @eval begin
                mutable struct MyCustomArray{T,N} <: AbstractArray{T,N}
                    data::Array{T,N}
                    sizehint_called::Bool
                end
                MyCustomArray(data::Array{T,N}) where {T,N} = MyCustomArray(data, false)
                Base.sizehint!(arr::MyCustomArray, args...) = (arr.sizehint_called = true; arr)
            end
            arr = QuantityArray(MyCustomArray([1.0, 2, 3]), dimension(u"m"));
            @test !ustrip(arr).sizehint_called
            sizehint!(arr, 5)
            @test ustrip(arr).sizehint_called
        end

        @testset "Reduce and vcat $Q" begin
            arr1 = QuantityArray([Q(us"hr")])
            arr2 = QuantityArray([Q(us"hr")])

            arr3 = vcat(arr1, arr2)
            @test typeof(arr3) <: QuantityArray{Float64,1,<:SymbolicDimensions,<:Q{Float64,<:SymbolicDimensions},<:Array}
            @test arr3 == QuantityArray(Q.([us"hr", us"hr"]))
            @test uexpand(arr3) == QuantityArray(Q.([u"hr", u"hr"]))

            arr4 = reduce(vcat, [arr1, arr2])
            @test typeof(arr4) <: QuantityArray{Float64,1,<:SymbolicDimensions,<:Q{Float64,<:SymbolicDimensions},<:Array}
            @test arr3 == arr4
            @test all(arr3 .== arr4)

            # Should work for incompatible units as well
            @test reduce(vcat, [[0.5us"hr"], [0.5us"m"]]) == [0.5u"hr", 0.5u"m"]
            @test typeof(reduce(vcat, [[0.5us"hr"], [0.5us"m"]])) <: Vector{<:Quantity{Float64,<:SymbolicDimensions}}
        end

        Q in (Quantity, RealQuantity) && @testset "Broadcast different arrays $Q" begin
            f(x, y, z, w) = x * y + z * w
            g(x, y, z, w) = f.(x, y, z, w)

            x = randn(32)
            y = QuantityArray(randn(32), u"km/s")
            z = rand(1:10, 32)
            w = Q{Float32}(u"m/s")
            @test typeof(g(x, y, z, w)) <: QuantityArray{Float64}

            y32 = QuantityArray(ustrip(y), dimension(y))
            @test typeof(y .* y32) <: QuantityArray{Float64}

            a = [randn() * u"km/s" for i=1:32]
            @test typeof(y .* a) <: QuantityArray
            @test typeof(a .* y) <: QuantityArray

            b = GenericQuantity(randn(Float32, 32), length=1, time=-1)
            @test typeof(b) <: GenericQuantity
            @test typeof(b .* b) <: Vector{<:GenericQuantity}
            @test typeof(a .* b) <: Vector{<:GenericQuantity}
            @test typeof(b .* a) <: Vector{<:GenericQuantity}
            @test typeof(y .* b) <: QuantityArray{Float64}
            @test typeof(b .* y) <: QuantityArray{Float64}
        end

        Q in (RealQuantity, Quantity) && @testset "Broadcast scalars $Q" begin
            for (x, qx) in ((0.5, 0.5u"s"), ([0.5, 0.2], GenericQuantity([0.5, 0.2], time=1)))
                @test size(qx) == size(x)
                @test length(qx) == length(x)
                @test axes(qx) == axes(x)
                @test iterate(qx)[1] == (iterate(x)[1] * u"s")
                @test ndims(qx) == ndims(x)
                @test Base.broadcastable(qx) == qx
                ustrip(qx) isa Real && @test qx[1] == qx
                @test keys(qx) == keys(x)
            end
        end

        @testset "Symbolic units $Q" begin
            z_ar = randn(32)
            z = QuantityArray(z_ar, Q(us"Constants.h * km/s"))
            z_expanded = QuantityArray(z_ar .* Q(u"Constants.h * km/s"))
            @test typeof(uexpand(z)) == typeof(z_expanded)
            @test all(uexpand(z) .≈ z_expanded)
            io = IOBuffer()
            Base.showarg(io, z, true)
            msg = String(take!(io))
            Q == Quantity && @test occursin(r"QuantityArray\(::Vector{Float64}, ::(DynamicQuantities\.)?Quantity{Float64, (DynamicQuantities\.)?SymbolicDimensions{(DynamicQuantities\.)?FixedRational{Int32, 25200}}}\)", msg)

            io = IOBuffer()
            Base.show(io, MIME"text/plain"(), typeof(z))
            msg2 = String(take!(io))
            @test msg2 == msg
        end

        Q == Quantity && @testset "Extra test coverage $Q" begin
            @test_throws ErrorException DynamicQuantities.materialize_first(())
            @test_throws "Unexpected broadcast" DynamicQuantities.materialize_first(())

            # Not sure how to test this otherwise, but method is supposed to be
            # required for the broadcasting interface
            x = [1u"km/s"]
            ref = Base.RefValue(x)
            @test DynamicQuantities.materialize_first(ref) === x[1]
        end
    end

    @testset "Unimplemented methods" begin
        qa = QuantityArray(randn(3), u"km/s")
        @test_throws ErrorException Base._similar_for(copy(qa), typeof(u"km/s"), qa, Base.SizeUnknown(), nothing)
        @test_throws ErrorException Base._similar_for(copy(qa), typeof(u"km/s"), qa, Base.HasLength(), 1)
        if hasmethod(Base._similar_for, Tuple{Array,Type,Any,Base.HasShape})
            @test_throws ErrorException Base._similar_for(copy(qa), typeof(u"km/s"), qa, Base.HasLength())
            @test_throws ErrorException Base._similar_for(copy(qa), typeof(u"km/s"), qa, Base.SizeUnknown())
            @test_throws ErrorException Base._similar_for(copy(qa), typeof(u"km/s"), qa, 1)
        end
    end
end

@testset "GenericQuantity" begin
    @testset "GenericQuantity construction" begin
        x = GenericQuantity(1.5)
        @test x isa GenericQuantity
        @test ustrip(x) == 1.5
        @test dimension(x) == Dimensions()

        x = GenericQuantity(big(1.5))
        @test typeof(x) <: GenericQuantity{BigFloat}

        x = GenericQuantity([1.5, 2.0], Dimensions{Rational{Int8}}; length=1)
        @test x isa GenericQuantity{Vector{Float64},Dimensions{Rational{Int8}}}
        @test supertype(typeof(x)) <: AbstractGenericQuantity{Vector{Float64},Dimensions{Rational{Int8}}}
        @test ustrip(x) == [1.5, 2.0]
        @test dimension(x) == Dimensions(length=1)

        x = GenericQuantity(randn(3,3))
        @test x isa GenericQuantity{Matrix{Float64}}
        @test size(x) == (3,3)

        x = GenericQuantity(1.5, length=1)
        y = GenericQuantity(2.0, time=1)
        @test_throws DimensionError x + y

        x = 0.5us"km/s"
        y = GenericQuantity(1.0)

    end

    @testset "GenericQuantity and Quantity promotion" begin
        x = GenericQuantity(1.5f0)
        y = Quantity(1.5, length=1)

        # *Always* promotes to GenericQuantity:
        @test typeof(x * y) <: GenericQuantity{Float64}

        x = GenericQuantity(rand(3))
        y = GenericQuantity(rand(3), length=1)

        @test typeof(x .* y) <: Vector{<:GenericQuantity{Float64}}

        x = 0.5us"km/s"
        @test GenericQuantity(x) isa GenericQuantity{Float64}

        x = GenericQuantity("abcd"; length=1)
        @test x isa GenericQuantity{String}
        y = "c"
        @test x * y == GenericQuantity("abcdc"; length=1)

        x = GenericQuantity([1.0, 2.0]; length=1)
        y = Quantity(3; mass=-1)
        @test x * y isa GenericQuantity{Vector{Float64}}
        @test x * y == GenericQuantity([3.0, 6.0]; length=1, mass=-1)

        x = [GenericQuantity(1.0; length=1), Quantity(2.0; length=1)]
        @test x isa Vector{<:GenericQuantity}
        ax = QuantityArray(x)
        @test ax isa QuantityArray{Float64,1,<:Dimensions,<:GenericQuantity{Float64}}

        # Test both symbolic/non-symbolic and Quantity/GenericQuantity:
        x = [0.5us"km/s", GenericQuantity(1.0f0; length=1)]
        @test x isa Vector{<:GenericQuantity{Float64,<:Dimensions}}
        @test ustrip(x[1]) == 500.0
        @test ustrip(x[2]) == 1.0

        x = [GenericQuantity([1.0, 2.0]), GenericQuantity([3f0, 4f0], Dimensions{Rational{Int}}, length=1)]
        @test x isa Vector{GenericQuantity{Vector{Float64},Dimensions{Rational{Int}}}}

        # Explicitly trigger Number conversion:
        x = Quantity(1u"nm")
        @test convert(Number, x) === x
        @test x isa Number
    end

    @testset "GenericQuantity broadcasting" begin
        x = QuantityArray([GenericQuantity(1.0f0i; length=1) for i=1:30])
        y = randn(30)
        z = QuantityArray([Quantity(1.0f0i; length=1) for i=1:30])

        @test x isa AbstractArray{<:GenericQuantity}
        @test x isa QuantityArray{Float32,1,<:Dimensions,<:GenericQuantity{Float32}}
        @test x .* y isa QuantityArray{Float64,1,<:Dimensions,<:GenericQuantity{Float64}}

        # TODO: Currently this converts to a `Vector` of `GenericQuantity`
        @test_skip x .* z isa QuantityArray{Float32,1,<:Dimensions,<:GenericQuantity{Float32}}
    end

    @testset "Array conversion" begin
        x = SArray{Tuple{3}}(randn(3))
        y = SArray{Tuple{3}}(randn(Float32, 3))
        qx = QuantityArray(x, Dimensions(Rational{Int}, length=1))
        qy = QuantityArray(y; length=1)

        @test typeof(convert(typeof(qx), qy)) == typeof(qx)
        @test convert(typeof(qx), qy)[1] isa Quantity{Float64}
        @test convert(typeof(qx), qy)[1] == convert(Quantity{Float64}, qy[1])
    end
end

function is_input_valid(f, x)
    try
        f(x)
    catch e
        e isa DomainError && return false
        rethrow(e)
    end
    return true
end

@testset "Assorted dimensionless functions" begin
    functions = (
        :sin, :cos, :tan, :sinh, :cosh, :tanh, :asin, :acos,
        :asinh, :acosh, :atanh, :sec, :csc, :cot, :asec, :acsc, :acot, :sech, :csch,
        :coth, :asech, :acsch, :acoth, :sinc, :cosc, :cosd, :cotd, :cscd, :secd,
        :sinpi, :cospi, :sind, :tand, :acosd, :acotd, :acscd, :asecd, :asind,
        :log, :log2, :log10, :log1p, :exp, :exp2, :exp10, :expm1, :frexp, :exponent,
        :atan, :atand
    )
    for Q in (RealQuantity, Quantity, GenericQuantity), D in (Dimensions, SymbolicDimensions), f in functions
        # Only test on valid domain
        valid_inputs = filter(
            x -> is_input_valid(eval(f), x),
            5rand(100) .- 2.5
        )
        for x in valid_inputs[1:3]
            qx_dimensionless = Q(x, D)
            qx_dimensions = convert(with_type_parameters(Q, Float64, D), Q(x, dimension(u"m/s")))
            @eval @test $f($qx_dimensionless) == $f($x)
            @eval @test_throws DimensionError $f($qx_dimensions)
            if f in (:atan, :atand)
                for y in valid_inputs[end-3:end]
                    qy_dimensionless = Q(y, D)
                    qy_dimensions = convert(with_type_parameters(Q, Float64, D), Q(y, dimension(u"m/s")))
                    @eval @test $f($y, $qx_dimensionless) == $f($y, $x)
                    @eval @test $f($qy_dimensionless, $x) == $f($y, $x)
                    @eval @test $f($qy_dimensionless, $qx_dimensionless) == $f($y, $x)
                    @eval @test $f($qy_dimensions, $qx_dimensions) == $f($y, $x)
                    @eval @test_throws DimensionError $f($qy_dimensions, $x)
                    @eval @test_throws DimensionError $f($y, $qx_dimensions)
                end
            end
        end
    end
    
    # Tests factorial (single integer input).
    for Q in (Quantity{Int64},), D in (Dimensions, SymbolicDimensions)
        for x in rand(1:10, 3)
            qx_dimensionless = Q(x, D)
            qx_dimensions = convert(with_type_parameters(Q, Int64, D), Q(x*u"m/s"))
            @eval @test $factorial($qx_dimensionless) == $factorial($x)
            @eval @test_throws DimensionError $factorial($qx_dimensions)
        end
    end

    # Tests binomial (two integer inputs).
    for Q in (Quantity{Int64},), D in (Dimensions, SymbolicDimensions)
        for x in rand(1:10, 3), y in rand(1:10, 3)
            qx_dimensionless = Q(x, D)
            qx_dimensions = convert(with_type_parameters(Q, Int64, D), Q(x*u"m/s"))
            qy_dimensionless = Q(y, D)
            qy_dimensions = convert(with_type_parameters(Q, Int64, D), Q(y*u"m/s"))
            @eval @test $binomial($y, $qx_dimensionless) == $binomial($y, $x)
            @eval @test $binomial($qy_dimensionless, $x) == $binomial($y, $x)
            @eval @test $binomial($qy_dimensionless, $qx_dimensionless) == $binomial($y, $x)
            @eval @test_throws DimensionError $binomial($qy_dimensions, $qx_dimensions) == $binomial($y, $x)
            @eval @test_throws DimensionError $binomial($qy_dimensions, $x)
            @eval @test_throws DimensionError $binomial($y, $qx_dimensions)
        end
    end

    s = record_show(DimensionError(u"km/s"), showerror)
    @test occursin("not dimensionless", s)
end

@testset "Assorted dimensionful functions" begin
    functions = (
        :float, :abs, :real, :imag, :conj, :adjoint, :unsigned,
        :nextfloat, :prevfloat, :identity, :transpose,
        :copysign, :flipsign, :modf,
        :floor, :trunc, :ceil, :significand,
        :ldexp, :round, :mod, :rem
    )
    for Q in (RealQuantity, Quantity, GenericQuantity), D in (Dimensions, SymbolicDimensions), f in functions
        T = f in (:abs, :real, :imag, :conj) ? ComplexF64 : Float64
        T <: Complex && Q == RealQuantity && continue
        if f == :modf  # Functions that return multiple outputs
            for x in 5rand(T, 3) .- 2.5
                qx_dimensions = convert(with_type_parameters(Q, T, D), Q(x, dimension(u"m/s")))
                num_outputs = 2
                for i=1:num_outputs
                    @eval @test $f($qx_dimensions)[$i] == $Q($f($x)[$i], $(dimension(u"m/s")))
                end
            end
        elseif f in (:copysign, :flipsign, :rem, :mod)  # Functions that need multiple inputs
            for x in 5rand(T, 3) .- 2.5
                for y in 5rand(T, 3) .- 2.5
                    # dim = convert(D, dimension(u"m/s"))
                    # qx_dimensions = Q(x, dim)
                    # qy_dimensions = Q(y, dim)
                    qx_dimensions = convert(with_type_parameters(Q, T, D), Q(x, dimension(u"m/s")))
                    qy_dimensions = convert(with_type_parameters(Q, T, D), Q(y, dimension(u"m/s")))
                    @eval @test $f($qx_dimensions, $qy_dimensions) == $Q($f($x, $y), dimension(u"m/s"))
                    if f in (:copysign, :flipsign)
                        # Also do test without dimensions
                        @eval @test $f($x, $qy_dimensions) == $f($x, $y)
                        @eval @test $f($qx_dimensions, $y) == $Q($f($x, $y), dimension(u"m/s"))
                    elseif f in (:rem, :mod)
                        # Also do test without dimensions (need dimensionless)
                        qx_dimensionless = Q(x, D)
                        qy_dimensionless = Q(y, D)
                        @eval @test $f($x, $qy_dimensionless) ≈ $Q($f($x, $y), $D)
                        @eval @test $f($qx_dimensionless, $y) ≈ $Q($f($x, $y), $D)
                        @eval @test_throws DimensionError $f($qx_dimensions, $y)
                        @eval @test_throws DimensionError $f($x, $qy_dimensions)
                        if f == :rem
                            # Can also do other rounding modes
                            for r in (:RoundFromZero, :RoundNearest, :RoundUp, :RoundDown)
                                @eval @test $f($qx_dimensions, $qy_dimensions, $r) ≈ $Q($f($x, $y, $r), dimension(u"m/s"))
                            end
                        end
                    end
                end
            end
        elseif f == :unsigned
            for x in 5rand(10:50, 3)
                qx_dimensions = convert(with_type_parameters(Q, typeof(x), D), Q(x, dimension(u"m/s")))
                @eval @test $f($qx_dimensions) == $Q($f($x), dimension(u"m/s"))
            end
        elseif f in (:round, :floor, :trunc, :ceil)
            for x in 5rand(T, 3) .- 2.5
                qx_dimensions = convert(with_type_parameters(Q, T, D), Q(x, dimension(u"m/s")))
                @eval @test $f($qx_dimensions) == $Q($f($x), dimension(u"m/s"))
                @eval @test $f(Int32, $qx_dimensions) == $Q($f(Int32, $x), dimension(u"m/s"))
            end
        elseif f == :ldexp
            for x in 5rand(T, 3) .- 2.5
                qx_dimensions = convert(with_type_parameters(Q, T, D), Q(x, dimension(u"m/s")))
                for i=1:3
                    @eval @test $f($qx_dimensions, $i) == $Q($f($x, $i), dimension(u"m/s"))
                end
            end
        else
            # Only test on valid domain
            valid_inputs = filter(
                x -> is_input_valid(eval(f), x),
                5rand(T, 100) .- 2.5
            )
            for x in valid_inputs[1:3]
                qx_dimensions = convert(with_type_parameters(Q, T, D), Q(x, dimension(u"m/s")))
                @eval @test $f($qx_dimensions) == $Q($f($x), dimension(u"m/s"))
            end
        end
    end
end

@testset "Assorted comparison functions" begin
    functions = (
        :(<=), :(<), :(>=), :(>), :isless, :isequal, :(==),
    )
    x = 5randn(10) .- 2.5
    y = 5randn(10) .- 2.5
    for Q in (RealQuantity, Quantity, GenericQuantity), D in (Dimensions, SymbolicDimensions), f in functions
        ground_truth = @eval $f.($x, $y)
        qx_dimensions = [convert(with_type_parameters(Q, Float64, D), Q(xi, dimension(u"m/s"))) for xi in x]
        qy_dimensions = [convert(with_type_parameters(Q, Float64, D), Q(yi, dimension(u"m/s"))) for yi in y]
        @eval @test all($f.($qx_dimensions, $qy_dimensions) .== $ground_truth)
        if f in (:isequal, :(==))
            # These include a dimension check in the result, rather than
            # throwing an error
            @eval @test !any($f.($qx_dimensions, $y))
            @eval @test !any($f.($x, $qy_dimensions))
        else
            @eval @test_throws DimensionError $f($qx_dimensions[1], $y[1])
            @eval @test_throws DimensionError $f($x[1], $qy_dimensions[1])
        end
        qx_dimensionless = [Q(xi, D) for xi in x]
        qy_dimensionless = [Q(yi, D) for yi in y]
        @eval @test all($f.($qx_dimensionless, $y) .== $ground_truth)
        @eval @test all($f.($x, $qy_dimensionless) .== $ground_truth)

        qx_real_dimensions = [convert(RealQuantity{Float64,D}, Quantity(xi, dimension(u"m/s"))) for xi in x]
        qy_real_dimensions = [convert(RealQuantity{Float64,D}, Quantity(yi, dimension(u"m/s"))) for yi in y]
        # Mixed quantity input
        @eval @test all($f.($qx_real_dimensions, $qy_dimensions) .== $ground_truth)
        @eval @test all($f.($qx_dimensions, $qy_real_dimensions) .== $ground_truth)
    end

    # Should be able to compare against `NoDims`:
    @test Quantity(1.0) >= 1.0
    @test !(Quantity(1.0) > 1.0)
end

@testset "Tests of `NoDims`" begin
    @test promote_type(NoDims{Int16}, NoDims{Int32}) === NoDims{Int32}

    # Prefer other types, always:
    @test promote_type(Dimensions{Int16}, NoDims{Int32}) === Dimensions{Int16}
    @test promote_type(MyDimensions{Int16}, NoDims{Int32}) === MyDimensions{Int16}

    # Always zero dimensions
    @test iszero(dimension(1.0))
    @test iszero(dimension([1.0, 2.0]))
    @test dimension(1.0) * u"1" == u"1"
    @test typeof(dimension(1.0) * u"1") === typeof(u"1")

    # Even when accessed:
    @test NoDims().km == 0
    @test NoDims().m != 1

    # Even weird user-defined dimensions:
    @test NoDims().cookie == 0

    # Always returns the same type:
    @test NoDims{Int32}().cookie isa Int32
    @test NoDims{Int32}().cookie == 0
end

@testset "Tests of SymbolicDimensionsSingleton" begin
    km = SymbolicUnits.km
    m = SymbolicUnits.m
    @test km isa Quantity{T,SymbolicDimensionsSingleton{R}} where {T,R}
    @test dimension(km) isa SymbolicDimensionsSingleton
    @test dimension(km) isa AbstractSymbolicDimensions

    @test dimension(km).km == 1
    @test dimension(km).m == 0
    @test_throws "is not available as a symbol" dimension(km).γ
    @test !iszero(dimension(km))
    @test inv(km) == us"km^-1"
    @test inv(km) == u"km^-1"

    @test !iszero(dimension(SymbolicConstants.c))
    @test SymbolicConstants.c isa Quantity{T,SymbolicDimensionsSingleton{R}} where {T,R}

    # Constructors
    @test SymbolicDimensionsSingleton(:cm) isa SymbolicDimensionsSingleton{DEFAULT_DIM_BASE_TYPE}
    @test constructorof(SymbolicDimensionsSingleton) === SymbolicDimensionsSingleton

    @test with_type_parameters(
            SymbolicDimensionsSingleton{Int64},
            Int32
        ) === SymbolicDimensionsSingleton{Int32}

    @test convert(
            SymbolicDimensions,
            SymbolicDimensionsSingleton{Int32}(:cm)
        ) isa SymbolicDimensions{Int32}

    @test copy(km) == km
    # Any operation should immediately convert it:
    @test km ^ -1 isa Quantity{T,DynamicQuantities.SymbolicDimensions{R}} where {T,R}

    # Test promotion explicitly for coverage:
    @test promote_type(
            SymbolicDimensionsSingleton{Int16},
            SymbolicDimensionsSingleton{Int32}
        ) === SymbolicDimensions{Int32}
    # ^ Note how we ALWAYS convert to SymbolicDimensions, even
    # if the types are the same.
    @test promote_type(
            SymbolicDimensionsSingleton{Int16},
            SymbolicDimensions{Int32}
        ) === SymbolicDimensions{Int32}
    @test promote_type(
            SymbolicDimensionsSingleton{Int64},
            Dimensions{Int16}
        ) === Dimensions{Int64}

    # Test map_dimensions explicitly for coverage:
    @test map_dimensions(-, dimension(km)).km == -1
    @test map_dimensions(-, dimension(km)) isa SymbolicDimensions
    @test map_dimensions(+, dimension(km), dimension(m)).km == 1
    @test map_dimensions(+, dimension(km), dimension(m)).m == 1
    @test map_dimensions(+, dimension(km), dimension(m)).cm == 0
    @test map_dimensions(+, dimension(km), SymbolicDimensions(dimension(m))).km == 1
    @test map_dimensions(+, dimension(km), SymbolicDimensions(dimension(m))).m == 1
    @test map_dimensions(+, dimension(km), SymbolicDimensions(dimension(m))).cm == 0
    @test map_dimensions(+, SymbolicDimensions(dimension(km)), dimension(m)).km == 1
    @test map_dimensions(+, SymbolicDimensions(dimension(km)), dimension(m)).m == 1
    @test map_dimensions(+, SymbolicDimensions(dimension(km)), dimension(m)).cm == 0

    # Note that we avoid converting to SymbolicDimensionsSingleton for uconvert:
    @test km |> uconvert(us"m") == 1000m
    @test km |> uconvert(us"m") isa Quantity{T,SymbolicDimensions{R}} where {T,R}
    @test [km, km] isa Vector{Quantity{T,SymbolicDimensionsSingleton{R}}} where {T,R}
    @test [km^2, km] isa Vector{Quantity{T,SymbolicDimensions{R}}} where {T,R}

    @test km |> uconvert(us"m") == km |> us"m"
    # No issue when converting to SymbolicDimensionsSingleton (gets
    # converted)
    @test uconvert(km, u"m") == 0.001km
    @test uconvert(km, u"m") isa Quantity{T,SymbolicDimensions{R}} where {T,R}

    # Symbolic dimensions retain symbols:
    @test QuantityArray([km, km]) |> uconvert(us"m") == [1000m, 1000m]
    @test QuantityArray([km, km]) |> uconvert(us"m") != [km, km]
end


@testset "Tests of AffineDimensions" begin
    °C  = ua"°C"
    °F  = ua"°F"
    mps = ua"m/s"

    @test aff_uparse("m/(s^2.5)") == ua"m/(s^2.5)"
    @test_throws ArgumentError aff_uparse("s[1]")
    @test_throws ArgumentError aff_uparse("pounds_per_hour")
    @test °C isa Quantity{T,AffineDimensions{R}} where {T,R}
    @test dimension(°C) isa AffineDimensions
    @test dimension(°C) isa AbstractAffineDimensions

    @test DynamicQuantities.basedim(dimension(°C)).temperature == 1
    @test DynamicQuantities.basedim(dimension(°C)).length == 0

    @test inv(mps) == us"s/m"
    @test inv(mps) == u"s/m"
    @test mps^2 == u"m^2/s^2"

    @test °C == ua"degC"
    @test °F == ua"degF"
    @test dimension(°C) == dimension(ua"degC")
    @test (°C - ua"degC") == 0.0u"K"

    # Constructors
    @test with_type_parameters(AffineDimensions, Float64) == AffineDimensions{Float64}
    @test constructorof(AffineDimensions) == AffineDimensions{DynamicQuantities.DEFAULT_DIM_BASE_TYPE}
    @test constructorof(AffineDimensions{Float64}) == AffineDimensions{Float64}
    @test Quantity(1.0, AffineDimensions(dimension(u"K"))) == u"K"
    @test AffineDimensions(scale=1, offset=0, basedim=dimension(u"K")) == AffineDimensions(basedim=dimension(u"K"))
    @test AffineDimensions(scale=1, offset=0, basedim=u"K") == AffineDimensions(basedim=dimension(ua"K"))
    @test AffineDimensions(scale=1.0, offset=273.15u"K", basedim=dimension(u"K")) == AffineDimensions(basedim=ua"°C")

    kelvin  = AffineDimensions(basedim=u"K")
    @test Quantity(1.0, kelvin) == u"K"

    rankine = AffineDimensions(scale=5/9, offset=0.0, basedim=dimension(u"K"))
    @test Quantity(1.0, rankine) == (5/9)u"K"

    fahrenheit = AffineDimensions(scale=1.0, offset=Quantity(459.67, rankine), basedim=rankine)
    @test Quantity(1.0, fahrenheit) ≈ °F

    celsius = AffineDimensions(scale=9/5, offset=Quantity(32.0, rankine), basedim=°F)
    @test Quantity(1.0, celsius) ≈ °C

    # Round-trip sanity checks
    @test -40°C ≈ -40°F
    @test Quantity(-40.0, celsius) ≈ Quantity(-40.0, fahrenheit)
    
    # Test promotion explicitly for coverage:
    @test promote_type(AffineDimensions{Int16}, AffineDimensions{Int32}) === AffineDimensions{Int32}
    @test promote_type(Dimensions{Int16}, AffineDimensions{Int32}) === Dimensions{Int32}
    @test promote_type(AffineDimensions{Int16}, Dimensions{Int32}) === Dimensions{Int32}
    @test promote_type(SymbolicDimensions{Int16}, AffineDimensions{Int32}) === Dimensions{Int32}
    @test promote_type(AffineDimensions{Int16}, SymbolicDimensions{Int32}) === Dimensions{Int32}

    # Type conversions
    @test convert(Quantity{Float64, AffineDimensions}, u"kg") isa Quantity{Float64, AffineDimensions{DynamicQuantities.DEFAULT_DIM_BASE_TYPE}}
    @test convert(Quantity{Float64, AffineDimensions{Float64}}, u"kg") isa Quantity{Float64, AffineDimensions{Float64}}
    @test convert(Quantity{Float64, Dimensions}, ua"kg") isa Quantity{Float64, Dimensions{DynamicQuantities.DEFAULT_DIM_BASE_TYPE}}
    
    # Test uncovered operations
    @test (2.0ua"m")^2 == (2.0u"m")^2
    @test dimension(ua"m")^Int32(2) == dimension(ua"m^2")
    @test 2.0u"m" + 2.0ua"m" === 4.0u"m"
    @test 2.0ua"m" + 2.0ua"m" === 4.0u"m"
    @test 2.0u"m" - 2.0ua"m" === 0.0u"m"
    @test 2.0ua"m" - 2.0ua"cm" === 1.98u"m"
    @test 5.0°C - 4.0°C === 1.0u"K"
    @test 2.0u"K" ≈ 2.0ua"K"
    @test 2.0ua"K" ≈ 2.0ua"K"
    @test 2.0ua"K" ≈ 2.0u"K"
    @test_throws AssertionError (2ua"°C")^2
    @test uexpand(2ua"°C") == 275.15u"K"

    # Test conversions
    @test °C |> us"K" isa Quantity{<:Real, <:SymbolicDimensions}
    @test 0°C |> us"K" == 273.15us"K"
    @test us"K" |> °C isa Quantity{<:Real, <:AffineDimensions}
    @test 0us"K" |> °C  == -273.15°C
    @test °C |> °F isa Quantity{<:Real, <:AffineDimensions}
    @test 0°C |> °F == 32°F

    @test QuantityArray([0,1]°C) |> uconvert(°F) isa QuantityArray{T, <:Any, AffineDimensions{R}} where {T,R}
    @test DynamicQuantities.affine_quantity(us"kPa") == u"kPa"

    # Test display against errors
    celsius = AffineDimensions(offset=273.15, basedim=u"K")
    psi = AffineDimensions(basedim=6.89476us"kPa")
    io = IOBuffer()
    @test isnothing(show(io, (dimension(°F), dimension(ua"K"), psi, celsius, fahrenheit)))

    # Test updating affine units
    @test_warn "unit °C already exists, skipping" DynamicQuantities.update_external_affine_unit(:°C, °C)
    @test_warn "unit °C already exists, skipping" DynamicQuantities.update_external_affine_unit(:°C, dimension(°C))
    @test_warn "unit °C already exists, skipping" DynamicQuantities.update_external_affine_unit(dimension(°C))
    @test_throws "Cannot register affine dimension if symbol is :nothing" DynamicQuantities.update_external_affine_unit(celsius)

end



@testset "Test div" begin
    for Q in (RealQuantity, Quantity, GenericQuantity)
        x = Q{Int}(10, length=1)
        y = Q{Int}(3, mass=-1)
        @test div(x, y) == Q{Int}(3, length=1, mass=1)
        @test div(x, 3) == Q{Int}(3, length=1)
        @test div(10, y) == Q{Int}(3, mass=1)
        @test div(x, y, RoundFromZero) == Q{Int}(4, length=1, mass=1)
        @test div(x, 3, RoundFromZero) == Q{Int}(4, length=1)
        @test div(10, y, RoundFromZero) == Q{Int}(4, mass=1)
    end
    # Also test mixed quantities:
    x = RealQuantity{Int}(10, length=1)
    y = Quantity{Int}(3, mass=-1)
    @test div(x, y) == Quantity{Int}(3, length=1, mass=1)
    @test typeof(div(x, y)) <: Quantity{Int}
    @test div(x, y, RoundFromZero) == Quantity{Int}(4, length=1, mass=1)
end

@testset "Exponentiation" begin
    for Q in (RealQuantity, Quantity, GenericQuantity)
        x = Q(2.0, length=1)

        # Quantity raised to dimensionless quantity
        p = Q(3)
        @test x^p == Q(8.0, length=3)

        x = Q(2.0, length=1)
        @test x^-2.0 == Q(0.25, length=-2)

        x = Q(2.0, time=1)
        @test x^0.0 ≈ Q(1.0, time=0)

        x = Q(0.5, mass=2)
        @test x^3.0 == Q(0.125, mass=6)

        # Q raised to negative float
        x = Q(4.0, current=2)
        @test x^-2.5 ≈ Q(4.0^-2.5, current=-5.0)

        # Q raised to positive float
        x = Q(3.0, amount=1)
        @test x^1.5 ≈ Q(5.196152422706632, amount=3//2)

        # Zero quantity edge cases
        x = Q(0.0, acceleration=1)
        @test x^1 == x
        @test x^0 == Q(1.0)

        # Complex number exponent
        if Q !== RealQuantity
            x = Q(0.7 + 0.2im)
            @test x^(1 - 3im) ≈ (0.7+0.2im)^(1 - 3im) * Q(1)
        end
    end

    x = RealQuantity(2.0)
    y = Quantity(2.0im)
    @test typeof(x^y) <: Quantity

    x = RealQuantity(2.0, length=1)
    y = Quantity(2.0im)
    @test_throws DimensionError x^y

    x = RealQuantity(2.0)
    y = Quantity(2.0im, mass=1)
    @test_throws DimensionError x^y
end

# `@testset` rewrites the test block with a `let...end`, resulting in an invalid
# local `const` (ref: src/units.jl:26). To avoid it, register units outside the
# test block.
map_count_before_registering = length(UNIT_MAPPING)
all_map_count_before_registering = length(ALL_MAPPING)
affine_count_before_registering = length(AFFINE_UNIT_MAPPING)

skipped_register_unit = false
#Registering Symbolic Units
if :MyV ∉ UNIT_SYMBOLS  # (In case we run this script twice)
    @eval @register_unit MyV u"V"
else
    skipped_register_unit = true
end
if :MySV ∉ UNIT_SYMBOLS
    @eval @register_unit MySV us"V"
end
if :MySV2 ∉ UNIT_SYMBOLS
    @eval @register_unit MySV2 us"km/h"
end
if :psi ∉ UNIT_SYMBOLS
    @eval @register_unit psi 6.89476us"kPa"
end

#Registering Affine Units
if :psig ∉ AFFINE_UNIT_SYMBOLS  #This example is in the documentation so it better work
    @eval @register_affine_unit psig AffineDimensions(offset=u"Constants.atm", basedim=u"psi")
else
    skipped_register_unit = true
end
if :My°C ∉ AFFINE_UNIT_SYMBOLS  # (In case we run this script twice)
    @eval @register_affine_unit My°C ua"°C"
end
if :My°C2 ∉ AFFINE_UNIT_SYMBOLS
    @eval @register_affine_unit My°C2 dimension(ua"°C")
end

@test_throws "Unit `m` is already defined as `1.0 m`" esc(_register_unit(:m, u"s"))
@test_throws "Unit `°C` is already defined as `1.0 °C`" esc(_register_affine_unit(:°C, ua"°C"))

# Constants as well:
@test_throws "Unit `Ryd` is already defined" esc(_register_unit(:Ryd, u"Constants.Ryd"))

@testset "Register Unit" begin
    MyV = u"MyV"
    MySV = u"MySV"
    MySV2 = u"MySV2"
    My°C = ua"My°C"
    My°C2 = ua"My°C2"
    psi = u"psi"
    psig = ua"psig"

    @test MyV === u"V"
    @test MyV == us"V"
    @test MySV == us"V"
    @test MySV2 == us"km/h"
    @test MySV == ua"V"
    @test MySV2 == ua"km/h"
    @test psi == ua"psi"
    @test psi == u"psi"
    @test psig == ua"psig"
    @test 0*psig == u"Constants.atm"
    @test My°C == ua"My°C"
    @test My°C == uexpand(ua"My°C")
    @test My°C2 == ua"My°C2"
    @test My°C2 == uexpand(ua"My°C2")

    if !skipped_register_unit
        @test length(UNIT_MAPPING) == map_count_before_registering + 4
        @test length(ALL_MAPPING) == all_map_count_before_registering + 4
        @test length(AFFINE_UNIT_MAPPING) == affine_count_before_registering + 7
    end

    for my_unit in (MySV, MyV)
        @test my_unit in UNIT_VALUES
        @test my_unit in ALL_VALUES
        @test my_unit in SYMBOLIC_UNIT_VALUES
        @test my_unit in AFFINE_UNIT_VALUES #Non-affine units should also be registered
    end
    
    for my_unit in (:MySV, :MyV)
        @test my_unit in UNIT_SYMBOLS
        @test my_unit in ALL_SYMBOLS
        @test my_unit in AFFINE_UNIT_SYMBOLS #Non-affine units should also be registered
    end

    for my_unit in (My°C, My°C2) #Affine units should only show up in the affine unit registry
        @test my_unit in AFFINE_UNIT_VALUES
    end

    for my_unit in (:My°C, :My°C2) #Affine units should only show up in the affine unit registry
        @test my_unit in AFFINE_UNIT_SYMBOLS
    end
end

push!(LOAD_PATH, joinpath(@__DIR__, "precompile_test"))

using ExternalUnitRegistration: MyWb
@testset "Type of External Unit" begin
    @test MyWb isa DEFAULT_QUANTITY_TYPE
    @test MyWb/u"m^2*kg*s^-2*A^-1" == 1.0
end

pop!(LOAD_PATH)
