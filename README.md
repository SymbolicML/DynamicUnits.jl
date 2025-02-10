<div align="center">

<img src="https://github.com/SymbolicML/DynamicQuantities.jl/assets/7593028/10b1e6b8-f1c5-43bc-97e6-4ddb4c175293" width=500>


[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://symbolicml.org/DynamicQuantities.jl/dev/)
[![Build Status](https://github.com/SymbolicML/DynamicQuantities.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SymbolicML/DynamicQuantities.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://coveralls.io/repos/github/SymbolicML/DynamicQuantities.jl/badge.svg?branch=main)](https://coveralls.io/github/SymbolicML/DynamicQuantities.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

</div>
  
DynamicQuantities defines a simple statically-typed `Quantity` type for Julia.
Physical dimensions are stored as a *value*, as opposed to a parametric type, as in [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).
This can greatly improve both runtime performance, by avoiding type instabilities, and startup time, as it avoids overspecializing methods.

- [Performance](#performance)
- [Usage](#usage)
  - [Constants](#constants)
  - [Symbolic Units](#symbolic-units)
  - [Arrays](#arrays)
  - [Unitful](#unitful)
- [Types](#types)
- [Vectors](#vectors)

## Performance

DynamicQuantities can greatly outperform Unitful
when the compiler cannot infer dimensions in a function:

```julia
julia> using BenchmarkTools, DynamicQuantities; import Unitful

julia> dyn_uni = 0.2u"m/s"
0.2 m s⁻¹

julia> unitful = convert(Unitful.Quantity, dyn_uni)
0.2 m s⁻¹

julia> f(x, i) = x ^ i * 0.3;

julia> @btime f($dyn_uni, 1);
  2.708 ns (0 allocations: 0 bytes)

julia> @btime f($unitful, 1);
  2.597 μs (30 allocations: 1.33 KiB)
```

**Note the μ and n: this is a 1000x speedup!**
Here, the DynamicQuantities quantity object allows the compiler to build a function that is type stable,
while the Unitful quantity object, which stores its dimensions in the type, requires type inference at runtime.

However, if the dimensions in your function *can* be inferred by the compiler,
then you can get better speeds with Unitful:

```julia
julia> g(x) = x ^ 2 * 0.3;

julia> @btime g($dyn_uni);
  1.791 ns (0 allocations: 0 bytes)

julia> @btime g($unitful);
  1.500 ns (0 allocations: 0 bytes)
```

While both of these are type stable,
because Unitful parametrizes the type on the dimensions, functions can specialize
to units and the compiler can optimize away units from the code.

## Usage

You can create a `Quantity` object 
by using the convenience macro `u"..."`:

```julia
julia> x = 0.3u"km/s"
300.0 m s⁻¹

julia> y = 42 * u"kg"
42.0 kg
```

or by importing explicitly:

```julia
julia> using DynamicQuantities: kPa

julia> room_temp = 100kPa
100000.0 m⁻¹ kg s⁻²
```

Note that `Units` is an exported submodule, so you can
also access this as `Units.kPa`. You may like to define

```julia
julia> const U = Units
```

so that you can simply write, say, `U.kPa` or `C.m_e`.

This supports a wide range of SI base and derived units, with common
prefixes.

You can also construct values explicitly with the `Quantity` type,
with a value and keyword arguments for the powers of the physical dimensions
(`mass`, `length`, `time`, `current`, `temperature`, `luminosity`, `amount`):

```julia
julia> x = Quantity(300.0, length=1, time=-1)
300.0 m s⁻¹
```

Elementary calculations with `+, -, *, /, ^, sqrt, cbrt, abs` are supported:

```julia
julia> x * y
12600.0 m kg s⁻¹

julia> x / y
7.142857142857143 m kg⁻¹ s⁻¹

julia> x ^ 3
2.7e7 m³ s⁻³

julia> x ^ -1
0.0033333333333333335 m⁻¹ s

julia> sqrt(x)
17.320508075688775 m¹ᐟ² s⁻¹ᐟ²

julia> x ^ 1.5
5196.152422706632 m³ᐟ² s⁻³ᐟ²
```

Each of these values has the same type, which means we don't need to perform type inference at runtime.

Furthermore, we can do dimensional analysis by detecting `DimensionError`:

```julia
julia> x + 3 * x
1.2 m¹ᐟ² kg

julia> x + y
ERROR: DimensionError: 0.3 m¹ᐟ² kg and 10.2 kg² s⁻² have incompatible dimensions
```

The dimensions of a `Quantity` can be accessed either with `dimension(quantity)` for the entire `Dimensions` object:

```julia
julia> dimension(x)
m¹ᐟ² kg
```

or with `umass`, `ulength`, etc., for the various dimensions:

```julia
julia> umass(x)
1//1

julia> ulength(x)
1//2
```

Finally, you can strip units with `ustrip`:
    
```julia
julia> ustrip(x)
0.2
```

### Constants

There are a variety of physical constants accessible
via the `Constants` submodule:

```julia
julia> Constants.c
2.99792458e8 m s⁻¹
```

which you may like to define as

```julia
julia> const C = Constants
```

These can also be used inside the `u"..."` macro:

```julia
julia> u"Constants.c * Hz"
2.99792458e8 m s⁻²
```

Similarly, you can just import each individual constant:

```julia
julia> using DynamicQuantities.Constants: h
```

For the full list, see the [docs](https://symbolicml.org/DynamicQuantities.jl/dev/constants/).


### Symbolic Units

You can also choose to not eagerly convert to SI base units,
instead leaving the units as the user had written them.
For example:

```julia
julia> q = 100us"cm * kPa"
100.0 cm kPa

julia> q^2
10000.0 cm² kPa²
```

You can convert to regular SI base units with
`uexpand`:

```julia
julia> uexpand(q^2)
1.0e6 kg² s⁻⁴
```

This also works with constants:

```julia
julia> x = us"Constants.c * Hz"
1.0 Hz c

julia> x^2
1.0 Hz² c²

julia> uexpand(x^2)
8.987551787368176e16 m² s⁻⁴
```

You can also convert a quantity in regular base SI units to symbolic units with the `|>` infix operator
```julia
julia> 5e-9u"m" |> us"nm"
5.0 nm
```

You can also convert between different symbolic units.
(Note that you can write this more explicitly
with `uconvert(us"nm", 5e-9u"m")`.)


Finally, you can also import these directly:

```julia
julia> using DynamicQuantities.SymbolicUnits: cm
```

or constants:

```julia
julia> using DynamicQuantities.SymbolicConstants: h
```

Note that `SymbolicUnits` and `SymbolicConstants` are exported,
so you can simply access these as `SymbolicUnits.cm` and `SymbolicConstants.h`,
respectively.

#### Custom Units

You can create custom units with the `@register_unit` macro:

```julia
julia> @register_unit OneFiveV 1.5u"V"
```

and then use it in calculations normally:

```julia
julia> x = us"OneFiveV"
1.0 OneFiveV

julia> x * 10u"A" |> us"W"
15.0 W

julia> 3us"V" |> us"OneFiveV"
2.0 OneFiveV
```

### Affine units
Units that have an offset (such as °C = K + 273.15) are an unfortunate fact of life: they are used extensively but often result in ambiguous mathematical operations (many other packages, such as Unitful.jl only support limited operations for affine dimensions). `AffineDimensions` seeks to extend DynamicQuantities.jl to reduce dependence on Unitful.jl, and enable handling/converting such units in a flexible, type-stable manner.

`AffineDimensions` are a generalization of `Dimensions` and `SymbolicDimensions`. While SymbolicDimensions essentially add a scale to Dimensions, AffineDimensions will add both a scale and an offset. Verious constructors can be used to construct `AffineDimensions` from other dimensions.
```
kelvin  = AffineDimensions(basedim=u"K") #Assumes a scale of 1 and offset 0
rankine = AffineDimensions(scale=5/9, offset=0.0, basedim=dimension(u"K")) #Rankine is a scaled version of Kelvin, offset is assumed to be of units 'basedim'
fahrenheit = AffineDimensions(scale=1.0, offset=Quantity(459.67, rankine), basedim=rankine) #Its best to make offset a `Quantity` to be explicit
celsius = AffineDimensions(scale=9/5, offset=Quantity(32.0, rankine), basedim=fahrenheit) #When AffineDimensiosn are used, offset starts with basedim's offset
```
#### Custom affine units
To access units from the affine unit registry, the string macro `ua"..."` can be used. This macro will always return quantities with AffineDimensions, even if a non-affine unit is called (it will simply have an offset of 0). Because AffineDimensions are a generalization of SymbolicDimensions, the affine unit registry will mirror the symbolic unit registry.
```
@register_unit psi 6.89476us"kPa"
u"psi"
>> 6894.76 m⁻¹ kg s⁻²
us"psi"
>> 1.0 psi
ua"psi"
>> 1.0 psi
```
However, strictly affine units cannot belong to the symbolic registry, so a different macro must be used on an AffineDimension (or quantity thereof)
```
@register_affine_unit psig AffineDimensions(offset=u"Constants.atm", basedim=u"psi") #Gauge pressure implies atmospheric offset
ua"psig"
>> 1.0 psig
us"psig"
>> ERROR: LoadError: ArgumentError: Symbol psig not found in `Units` or `Constants`.
```
Affine unit parsing can also be done outside of a macro using `aff_uparse(str::AbstractString)`
```
aff_uparse("°C")
>> 1.0 °C
```
#### Operations on affine quantities
In Unitful.jl, multiplication of affine quantities is not supported for affine dimensions:
```
using Unitful
u"R"*0u"°C"
>> ERROR: AffineError: an invalid operation was attempted with affine units: °C
```
This behaviour is mimicked in DynamicQuantities:
```
using DynamicQuantities
u"Constants.R"*(0ua"°C")
>> AssertionError: AffineDimensions °C has a non-zero offset, implicit conversion is not allowed due to ambiguity. Use uexpand(x) to explicitly convert 
```
In general, it's best to treat quantities with AffineDimensions as placeholders and use `uexpand(q)` or `uconvert(units, q)` as soon as possible. The main objective of AffineDimesnions is to provide you with convenient, type-stable tools to do this conversion before applying mathematical operations.


### Arrays

For working with an array of quantities that have the same dimensions,
you can use a `QuantityArray`:

```julia
julia> ar = QuantityArray(rand(3), u"m/s")
3-element QuantityArray(::Vector{Float64}, ::Quantity{Float64, Dimensions{FixedRational{Int32, 25200}}}):
 0.2729202669351497 m s⁻¹
 0.992546340360901 m s⁻¹
 0.16863543422972482 m s⁻¹
```

This `QuantityArray` is a subtype `<:AbstractArray{Quantity{Float64,Dimensions{...}},1}`,
meaning that indexing a specific element will return a `Quantity`:

```julia
julia> ar[2]
0.992546340360901 m s⁻¹

julia> ar[2] *= 2
1.985092680721802 m s⁻¹

julia> ar[2] += 0.5u"m/s"
2.485092680721802 m s⁻¹
```

This also has a custom broadcasting interface which
allows the compiler to avoid redundant dimension calculations,
relative to if you had simply used an array of quantities:

```julia
julia> f(v) = v^2 * 1.5;

julia> @btime $f.(xa) setup=(xa = randn(100000) .* u"km/s");
  109.500 μs (2 allocations: 3.81 MiB)

julia> @btime $f.(qa) setup=(xa = randn(100000) .* u"km/s"; qa = QuantityArray(xa));
  50.917 μs (3 allocations: 781.34 KiB)
```

So we can see the `QuantityArray` version saves on both time and memory.

### Unitful

DynamicQuantities allows you to convert back and forth from Unitful.jl:

```julia
julia> using Unitful: Unitful, @u_str; import DynamicQuantities

julia> x = 0.5u"km/s"
0.5 km s⁻¹

julia> y = convert(DynamicQuantities.Quantity, x)
500.0 m s⁻¹

julia> y2 = y^2 * 0.3
75000.0 m² s⁻²

julia> x2 = convert(Unitful.Quantity, y2)
75000.0 m² s⁻²

julia> x^2*0.3 == x2
true
```

## Types

Both a `Quantity`'s values and dimensions are of arbitrary type. The default
`Dimensions` (for the `u"..."` macro) performs exponent tracking for SI units,
and `SymbolicDimensions` (for the `us"..."` macro) performs exponent tracking
for all known unit and constant symbols, using a sparse array.

You can create custom spaces dimension spaces by simply creating
a Julia struct subtyped to `AbstractDimensions`:

```julia
julia> struct CookiesAndMilk{R} <: AbstractDimensions{R}
           cookies::R
           milk::R
       end

julia> cookie_rate = Quantity(0.9, CookiesAndMilk(cookies=1, milk=-1))
0.9 cookies milk⁻¹

julia> total_milk = Quantity(103, CookiesAndMilk(milk=1))
103 milk

julia> total_cookies = cookie_rate * total_milk
92.7 cookies
```

Exponents are tracked by default with the type `R = FixedRational{Int32,C}`,
which represents rational numbers with a fixed denominator `C`.
This is much faster than `Rational`.

```julia
julia> typeof(0.5u"kg")
Quantity{Float64, Dimensions{FixedRational{Int32, 25200}}}
```

You can change the type of the value field by initializing with a value
explicitly of the desired type.

```julia
julia> typeof(Quantity(Float16(0.5), mass=1, length=1))
Quantity{Float16, Dimensions{FixedRational{Int32, 25200}}}
```

or by conversion:

```julia
julia> typeof(convert(Quantity{Float16}, 0.5u"m/s"))
Quantity{Float16, Dimensions{FixedRational{Int32, 25200}}}
```

For many applications, `FixedRational{Int8,6}` will suffice,
and can be faster as it means the entire `Dimensions`
struct will fit into 64 bits.
You can change the type of the dimensions field by passing
the type you wish to use as the second argument to `Quantity`:

```julia
julia> using DynamicQuantities

julia> R8 = Dimensions{FixedRational{Int8,6}};

julia> R32 = Dimensions{FixedRational{Int32,2^4 * 3^2 * 5^2 * 7}};  # Default

julia> q8 = [Quantity{Float64,R8}(randn(), length=rand(-2:2)) for i in 1:1000];

julia> q32 = [Quantity{Float64,R32}(randn(), length=rand(-2:2)) for i in 1:1000];

julia> f(x) = @. x ^ 2 * 0.5;

julia> @btime f($q8);
  1.433 μs (3 allocations: 15.77 KiB)

julia> @btime f($q32);
  1.883 μs (4 allocations: 39.12 KiB)
```
