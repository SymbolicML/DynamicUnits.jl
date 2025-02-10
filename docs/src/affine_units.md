# Affine Dimensions
Units that have an offset (such as °C = K + 273.15) are an unfortunate fact of life: they are used extensively but often result in ambiguous mathematical operations (many other packages, such as Unitful.jl only support limited operations for affine dimensions). `AffineDimensions` seeks to extend DynamicQuantities.jl to reduce dependence on Unitful.jl, and enable handling/converting such units in a flexible, type-stable manner.

`AffineDimensions` are a generalization of `Dimensions` and `SymbolicDimensions`. While SymbolicDimensions essentially add a scale to Dimensions, AffineDimensions will add both a scale and an offset. Verious constructors can be used to construct `AffineDimensions` from other dimensions.
```
kelvin  = AffineDimensions(basedim=u"K") #Assumes a scale of 1 and offset 0
rankine = AffineDimensions(scale=5/9, offset=0.0, basedim=dimension(u"K")) #Rankine is a scaled version of Kelvin, offset is assumed to be of units 'basedim'
fahrenheit = AffineDimensions(scale=1.0, offset=Quantity(459.67, rankine), basedim=rankine) #Its best to make offset a `Quantity` to be explicit
celsius = AffineDimensions(scale=9/5, offset=Quantity(32.0, rankine), basedim=fahrenheit) #When AffineDimensiosn are used, offset starts with basedim's offset
```
## Registration and parsing
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
```@docs
@ua_str
@register_affine_unit
aff_uparse
```
## Operations
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