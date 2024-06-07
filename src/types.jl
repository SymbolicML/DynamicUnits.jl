using Tricks: static_fieldnames

const DEFAULT_DIM_BASE_TYPE = FixedRational{DEFAULT_NUMERATOR_TYPE,DEFAULT_DENOM}
const DEFAULT_VALUE_TYPE = Float64

"""
    AbstractDimensions{R}

An abstract type for dimension types. `R` is the type of the exponents of the dimensions,
and by default is set to `DynamicQuantities.DEFAULT_DIM_BASE_TYPE`.
AbstractDimensions are used to store the dimensions of `UnionAbstractQuantity` objects.
Together these enable many operators in Base to manipulate dimensions.
This type has generic constructors for creating dimension objects, so user-defined
dimension types can be created by simply subtyping `AbstractDimensions`, without
the need to define many other functions.

The key function that one could wish to overload is
`DynamicQuantities.dimension_name(::AbstractDimensions, k::Symbol)` for mapping from a field name
to a base unit (e.g., `length` by default maps to `m`). You may also need to overload
`constructorof(::Type{T})` in case of non-standard construction.
"""
abstract type AbstractDimensions{R} end

"""
    AbstractQuantity{T,D} <: Number

An abstract type for quantities. `T` is the type of the value of the quantity,
which should be `<:Number`.
`D` is the type of the dimensions of the quantity. By default, `D` is set to
`DynamicQuantities.DEFAULT_DIM_TYPE`. `T` is inferred from the value in a calculation,
but in other cases is defaulted to `DynamicQuantities.DEFAULT_VALUE_TYPE`.
It is assumed that the value is stored in the `:value` field, and the dimensions
object is stored in the `:dimensions` field. These fields can be accessed with
`ustrip` and `dimension`, respectively. Many operators in `Base` are defined on
`AbstractQuantity` objects, including `+, -, *, /, ^, sqrt, cbrt, abs`.

See also `AbstractGenericQuantity` for creating quantities subtyped to `Any`.

**Note**: In general, you should probably
specialize on `UnionAbstractQuantity` which is
the union of both `AbstractQuantity` and `AbstractGenericQuantity`,
_as well as any other future abstract quantity types_,
"""
abstract type AbstractQuantity{T,D} <: Number end

"""
    AbstractGenericQuantity{T,D} <: Any

This has the same behavior as `AbstractQuantity` but is subtyped to `Any` rather
than `Number`.

**Note**: In general, you should probably
specialize on `UnionAbstractQuantity` which is
the union of both `AbstractQuantity` and `AbstractGenericQuantity`,
_as well as any other future abstract quantity types_,
"""
abstract type AbstractGenericQuantity{T,D} end

"""
    AbstractRealQuantity{T,D} <: Real

This has the same behavior as `AbstractQuantity` but is subtyped to `Real` rather
than `Number`.
"""
abstract type AbstractRealQuantity{T,D} <: Real end

"""
    UnionAbstractQuantity{T,D}

This is a union of both `AbstractQuantity{T,D}` and `AbstractGenericQuantity{T,D}`.
It is used throughout the library to declare methods which can take both types.
You should generally specialize on this type, rather than its constituents,
as it will also include future abstract quantity types.
"""
const UnionAbstractQuantity{T,D} = Union{AbstractQuantity{T,D},AbstractGenericQuantity{T,D},AbstractRealQuantity{T,D}}

"""
    Dimensions{R<:Real} <: AbstractDimensions{R}

A type representing the dimensions of a quantity, with each
field giving the power of the corresponding dimension. For
example, the dimensions of velocity are `Dimensions(length=1, time=-1)`.
Each of the 7 dimensions are stored using the type `R`,
which is by default a rational number.

# Fields

- `length`: length dimension (i.e., meters^(length))
- `mass`: mass dimension (i.e., kg^(mass))
- `time`: time dimension (i.e., s^(time))
- `current`: current dimension (i.e., A^(current))
- `temperature`: temperature dimension (i.e., K^(temperature))
- `luminosity`: luminosity dimension (i.e., cd^(luminosity))
- `amount`: amount dimension (i.e., mol^(amount))

# Constructors

- `Dimensions(args...)`: Pass all the dimensions as arguments.
- `Dimensions(; kws...)`: Pass a subset of dimensions as keyword arguments. `R` is set to `DEFAULT_DIM_BASE_TYPE`.
- `Dimensions(::Type{R}; kws...)` or `Dimensions{R}(; kws...)`: Pass a subset of dimensions as keyword arguments, with the output type set to `Dimensions{R}`.
- `Dimensions{R}()`: Create a dimensionless object typed as `Dimensions{R}`.
- `Dimensions{R}(d::Dimensions)`: Copy the dimensions from another `Dimensions` object, with the output type set to `Dimensions{R}`.
"""
struct Dimensions{R<:Real} <: AbstractDimensions{R}
    length::R
    mass::R
    time::R
    current::R
    temperature::R
    luminosity::R
    amount::R
end

(::Type{D})(::Type{R}; kws...) where {R,D<:AbstractDimensions} = with_type_parameters(D, R)((tryrationalize(R, get(kws, k, zero(R))) for k in dimension_names(D))...)
(::Type{D})(; kws...) where {R,D<:AbstractDimensions{R}} = constructorof(D)(R; kws...)
(::Type{D})(; kws...) where {D<:AbstractDimensions} = D(DEFAULT_DIM_BASE_TYPE; kws...)
function (::Type{D})(d::D2) where {R,D<:AbstractDimensions{R},D2<:AbstractDimensions}
    dimension_names_equal(D, D2) ||
        error("Cannot create a dimensions of `$(D)` from `$(D2)`. Please write a custom method for construction.")
    D((getproperty(d, k) for k in dimension_names(D))...)
end

const DEFAULT_DIM_TYPE = Dimensions{DEFAULT_DIM_BASE_TYPE}

"""
    NoDims{R}

A type representing the dimensions of a non-quantity.

For any `getproperty` call on this type, the result is `zero(R)`.
"""
struct NoDims{R<:Real} <: AbstractDimensions{R}
end

Base.getproperty(::NoDims{R}, ::Symbol) where {R} = zero(R)

const DEFAULT_DIMENSIONLESS_TYPE = NoDims{DEFAULT_DIM_BASE_TYPE}

"""
    Quantity{T<:Number,D<:AbstractDimensions} <: AbstractQuantity{T,D} <: Number

Physical quantity with value `value` of type `T` and dimensions `dimensions` of type `D`.
For example, the velocity of an object with mass 1 kg and velocity
2 m/s is `Quantity(2, mass=1, length=1, time=-1)`.
You should access these fields with `ustrip(q)`, and `dimension(q)`.
You can access specific dimensions with `ulength(q)`, `umass(q)`, `utime(q)`,
`ucurrent(q)`, `utemperature(q)`, `uluminosity(q)`, and `uamount(q)`.

Severals operators in `Base` are extended to work with `Quantity` objects,
including `*`, `+`, `-`, `/`, `abs`, `^`, `sqrt`, and `cbrt`, which manipulate
dimensions according to the operation.

# Fields

- `value::T`: value of the quantity of some type `T`. Access with `ustrip(::Quantity)`
- `dimensions::D`: dimensions of the quantity. Access with `dimension(::Quantity)`

# Constructors

- `Quantity(x; kws...)`: Construct a quantity with value `x` and dimensions given by the keyword arguments. The value
   type is inferred from `x`. `R` is set to `DEFAULT_DIM_TYPE`.
- `Quantity(x, ::Type{D}; kws...)`: Construct a quantity with value `x` with dimensions given by the keyword arguments,
   and the dimensions type set to `D`.
- `Quantity(x, d::D)`: Construct a quantity with value `x` and dimensions `d` of type `D`.
- `Quantity{T}(...)`: As above, but converting the value to type `T`. You may also pass a `Quantity` as input.
- `Quantity{T,D}(...)`: As above, but converting the value to type `T` and dimensions to `D`. You may also pass a
  `Quantity` as input.
"""
struct Quantity{T<:Number,D<:AbstractDimensions} <: AbstractQuantity{T,D}
    value::T
    dimensions::D

    Quantity(x::_T, dimensions::_D) where {_T,_D<:AbstractDimensions} = new{_T,_D}(x, dimensions)
end

"""
    GenericQuantity{T<:Any,D<:AbstractDimensions} <: AbstractGenericQuantity{T,D} <: Any

This has the same behavior as `Quantity` but is subtyped to `AbstractGenericQuantity <: Any`
rather than `AbstractQuantity <: Number`.
"""
struct GenericQuantity{T,D<:AbstractDimensions} <: AbstractGenericQuantity{T,D}
    value::T
    dimensions::D

    GenericQuantity(x::_T, dimensions::_D) where {_T,_D<:AbstractDimensions} = new{_T,_D}(x, dimensions)
end

"""
    RealQuantity{T<:Real,D<:AbstractDimensions} <: AbstractRealQuantity{T,D} <: Real

This has the same behavior as `Quantity` but is subtyped to `AbstractRealQuantity <: Real`.
"""
struct RealQuantity{T<:Real,D<:AbstractDimensions} <: AbstractRealQuantity{T,D}
    value::T
    dimensions::D

    RealQuantity(x::_T, dimensions::_D) where {_T,_D<:AbstractDimensions} = new{_T,_D}(x, dimensions)
end

"""
    ABSTRACT_QUANTITY_TYPES

A constant tuple of the existing abstract quantity types,
each as a tuple with (1) the abstract type,
(2) the base type, and (3) the default exported concrete type.
"""
const ABSTRACT_QUANTITY_TYPES = (
    (AbstractQuantity, Number, Quantity),
    (AbstractGenericQuantity, Any, GenericQuantity),
    (AbstractRealQuantity, Real, RealQuantity)
)


for (type, base_type, _) in ABSTRACT_QUANTITY_TYPES
    @eval begin
        (::Type{Q})(x::T, ::Type{D}; kws...) where {D<:AbstractDimensions,T<:$base_type,T2,Q<:$type{T2}} = constructorof(Q)(convert(T2, x), D(; kws...))
        (::Type{Q})(x::$base_type, ::Type{D}; kws...) where {D<:AbstractDimensions,Q<:$type} = constructorof(Q)(x, D(; kws...))
        (::Type{Q})(x::T; kws...) where {T<:$base_type,T2,Q<:$type{T2}} = constructorof(Q)(convert(T2, x), dim_type(Q)(; kws...))
        (::Type{Q})(x::$base_type; kws...) where {Q<:$type} = constructorof(Q)(x, dim_type(Q)(; kws...))
    end
    for (type2, _, _) in ABSTRACT_QUANTITY_TYPES
        @eval begin
            (::Type{Q})(q::$type2) where {T,D<:AbstractDimensions,Q<:$type{T,D}} = constructorof(Q)(convert(T, ustrip(q)), convert(D, dimension(q)))
            (::Type{Q})(q::$type2) where {T,Q<:$type{T}} = constructorof(Q)(convert(T, ustrip(q)), dimension(q))
            (::Type{Q})(q::$type2) where {Q<:$type} = constructorof(Q)(ustrip(q), dimension(q))
        end
    end
end

const DEFAULT_QUANTITY_TYPE = Quantity{DEFAULT_VALUE_TYPE, DEFAULT_DIM_TYPE}

@inline function new_dimensions(::Type{D}, dims...) where {D<:AbstractDimensions}
    return constructorof(D)(dims...)
end
@inline function new_quantity(::Type{Q}, val, dims) where {Q<:UnionAbstractQuantity}
    Qout = promote_quantity_on_value(Q, typeof(val))
    return constructorof(Qout)(val, dims)
end

Base.eltype(::Union{Type{<:AbstractDimensions{R}},AbstractDimensions{R}}) where {R} = R
dim_type(::Type{Q}) where {T,D<:AbstractDimensions,Q<:UnionAbstractQuantity{T,D}} = D
dim_type(::Type{<:UnionAbstractQuantity}) = DEFAULT_DIM_TYPE

"""
    constructorof(::Type{<:AbstractDimensions})
    constructorof(::Type{<:UnionAbstractQuantity})

Return the constructor of the given type. This is used to create new objects
of the same type as the input. Overload a method for a new type, especially
if you need custom behavior.
"""
constructorof(::Type{<:Dimensions}) = Dimensions
constructorof(::Type{<:Quantity}) = Quantity
constructorof(::Type{<:GenericQuantity}) = GenericQuantity
constructorof(::Type{<:RealQuantity}) = RealQuantity

"""
    with_type_parameters(::Type{<:AbstractDimensions}, ::Type{R})
    with_type_parameters(::Type{<:UnionAbstractQuantity}, ::Type{T}, ::Type{D})

Return the type with the given type parameters instead of the ones in the input type.
This is used to get `Dimensions{R}` from input `(Dimensions{R1}, R)`, for example.
Overload a method for a new type, especially if you need custom behavior.
"""
function with_type_parameters(::Type{<:Dimensions}, ::Type{R}) where {R}
    return Dimensions{R}
end
function with_type_parameters(::Type{<:Quantity}, ::Type{T}, ::Type{D}) where {T,D}
    return Quantity{T,D}
end
function with_type_parameters(::Type{<:GenericQuantity}, ::Type{T}, ::Type{D}) where {T,D}
    return GenericQuantity{T,D}
end
function with_type_parameters(::Type{<:RealQuantity}, ::Type{T}, ::Type{D}) where {T,D}
    return RealQuantity{T,D}
end

# The following functions should be overloaded for special types
function constructorof(::Type{T}) where {T<:Union{UnionAbstractQuantity,AbstractDimensions}}
    return Base.typename(T).wrapper
end
function with_type_parameters(::Type{D}, ::Type{R}) where {D<:AbstractDimensions,R}
    return constructorof(D){R}
end
function with_type_parameters(::Type{Q}, ::Type{T}, ::Type{D}) where {Q<:UnionAbstractQuantity,T,D}
    return constructorof(Q){T,D}
end

"""
    dimension_names(::Type{<:AbstractDimensions})

Return a tuple of symbols with the names of the dimensions of the given type.
This should be static so that it can be hardcoded during compilation.
The default is to use `fieldnames`, but you can overload this for custom behavior.
"""
@inline function dimension_names(::Type{D}) where {D<:AbstractDimensions}
    return static_fieldnames(D)
end

struct DimensionError{Q1,Q2} <: Exception
    q1::Q1
    q2::Q2

    DimensionError(q1, q2) = new{typeof(q1),typeof(q2)}(q1, q2)
    DimensionError(q1) = DimensionError(q1, nothing)
end
