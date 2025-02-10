
const AbstractQuantityOrArray{T,D} = Union{UnionAbstractQuantity{T,D}, QuantityArray{T,<:Any,D}}

abstract type AbstractAffineDimensions{R} <: AbstractDimensions{R} end

const AffineOrSymbolicDimensions{R} = Union{AbstractAffineDimensions{R}, AbstractSymbolicDimensions{R}}

@kwdef struct AffineDimensions{R} <: AbstractAffineDimensions{R}
    scale::Float64 = 1.0
    offset::Float64 = 0.0
    basedim::Dimensions{R}
    symbol::Symbol = :nothing
end

#Inferring the type parameter R ========================================================================================================================
AffineDimensions(s, o, dims::AbstractDimensions{R}, sym::Symbol=:nothing) where {R} = AffineDimensions{R}(s, o, dims, sym)
AffineDimensions(s, o, q::UnionAbstractQuantity{<:Any,<:AbstractDimensions{R}}, sym::Symbol=:nothing) where {R} = AffineDimensions{R}(s, o, q, sym)
AffineDimensions(d::Dimensions{R}) where R = AffineDimensions{R}(scale=1.0, offset=0.0, basedim=d, symbol=:nothing)

#Affine dimensions from other affine dimensions =========================================================================================================
function AffineDimensions{R}(s::Real, o::Real, dims::AbstractAffineDimensions, sym::Symbol=:nothing) where {R}
    new_s = s*scale(dims)
    new_o = offset(dims) + o*scale(dims) #Scale of o is assumed to be scale of base dimensions
    return AffineDimensions{R}(new_s, new_o, basedim(dims), sym)
end


function AffineDimensions{R}(s::Real, o::UnionAbstractQuantity, dims::AbstractAffineDimensions, sym::Symbol=:nothing) where {R}
    new_s = s*scale(dims)
    new_o = offset(dims) + ustrip(siunits(o)) #Offset is always in SI units
    return AffineDimensions{R}(new_s, new_o, basedim(dims), sym)
end

function AffineDimensions{R}(s::Real, o::UnionAbstractQuantity, dims::Dimensions, sym::Symbol=:nothing) where {R}
    return AffineDimensions{R}(s, ustrip(siunits(o)), dims, sym)
end

#Affine dimensions from quantities =========================================================================================================================
function AffineDimensions{R}(s::Real, o::UnionAbstractQuantity, q::UnionAbstractQuantity, sym::Symbol=:nothing) where R
    q0   = siunits(0*q) #Origin point in SI units
    oΔ   = siunits(o) - siunits(0*o) #Offset is a difference in affine units
    dimension(q0) == dimension(oΔ) || throw(DimensionError(o, q)) #Check the units and give an informative error
    
    #Obtain SI units of the scale and offset
    o_si = oΔ + q0 #Total offset is origin plus the offset
    q_si = siunits(q) - q0 #The scaling quantity must remove the origin
    
    #Call the SI quantity constructor
    return AffineDimensions{R}(s, o_si, q_si, sym)
end

#Base case when everyting is convrted to si units (offset is assumed to be in SI units)
function AffineDimensions{R}(s::Real, o::UnionAbstractQuantity{<:Any,<:Dimensions}, q::UnionAbstractQuantity{<:Any,<:Dimensions}, sym::Symbol=:nothing) where R
    dimension(o) == dimension(q) || throw(DimensionError(o, q))
    o_val = ustrip(o)
    q_val = ustrip(q)
    return AffineDimensions{R}(s*q_val, o_val, dimension(q), sym)
end

#If a quantity is used only for the dimension, the offset is assumed to be in the same scale as the quantity 
function AffineDimensions{R}(s::Real, o::Real, q::Q, sym::Symbol=:nothing) where {R, Q<:UnionAbstractQuantity}
    return AffineDimensions{R}(s, o*q, q, sym)
end


scale(d::AffineDimensions)  = d.scale
offset(d::AffineDimensions) = d.offset
basedim(d::AffineDimensions) = d.basedim

with_type_parameters(::Type{<:AffineDimensions}, ::Type{R}) where {R} = AffineDimensions{R}
constructorof(::Type{AffineDimensions}) = AffineDimensions{DEFAULT_DIM_BASE_TYPE}
constructorof(::Type{AffineDimensions{R}}) where R = AffineDimensions{R}

function Base.show(io::IO, d::AbstractAffineDimensions)
    addsign = ifelse(offset(d)<0, "-" , "+")

    if d.symbol != :nothing
        print(io, d.symbol)
    elseif isone(scale(d)) & iszero(offset(d))
        print(io, basedim(d))
    elseif iszero(offset(d))
        print(io, "(", scale(d), " ", basedim(d),")")
    elseif isone(scale(d))
        print(io, "(", addsign, abs(offset(d)), basedim(d), ")")
    else
        print(io, "(", scale(d), addsign, abs(offset(d)), " ", basedim(d),")")
    end
end

assert_no_offset(d::AffineDimensions) = iszero(offset(d)) || throw(AssertionError("AffineDimensions $(d) has a non-zero offset, implicit conversion is not allowed due to ambiguity. Use uexpand(x) to explicitly convert"))
siunits(q::UnionAbstractQuantity{<:Any, <:Dimensions}) = q
siunits(q::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}) = uexpand(q)
function siunits(q::Q) where {T, R, D<:AbstractAffineDimensions{R}, Q<:UnionAbstractQuantity{T,D}}
    return force_convert(with_type_parameters(Q, T, Dimensions{R}), q)
end
siunits(q::QuantityArray) = siunits.(q)


"""
uexpand(q::Q) where {T, R, D<:AbstractAffineDimensions{R}, Q<:UnionAbstractQuantity{T,D}}

Expand the affine units in a quantity to their base SI form. In other words, this converts a quantity with AbstractAffineDimensions   
to one with Dimensions. The opposite of this function is uconvert, for converting to specific symbolic units, or, e.g.,
convert(Quantity{<:Any,<:AbstractSymbolicDimensions}, q), for assuming SI units as the output symbols.
"""
uexpand(q::UnionAbstractQuantity{<:Any, <:AbstractAffineDimensions}) = siunits(q)


"""
affine_quantity(q::UnionAbstractQuantity)

Converts a quantity to its nearest affine quantity representation (with scale=1.0 and offset=0.0)
"""
function affine_quantity(q::Q) where {T, R, D<:AbstractDimensions{R}, Q<:UnionAbstractQuantity{T,D}}
    q_si  = siunits(q)
    dims  = AffineDimensions{R}(scale=1.0, offset=0.0, basedim=dimension(q_si))
    q_val = convert(T, ustrip(q_si))
    return constructorof(Q)(q_val, dims)
end

"""
affine_unit(q::UnionAbstractQuantity)

Converts a quantity to its nearest affine unit (with scale=ustrip(q) and offset=0.0)
"""
function affine_unit(q::Q) where {T, R, D<:AbstractDimensions{R}, Q<:UnionAbstractQuantity{T,D}}
    q_si  = siunits(q)
    dims  = AffineDimensions{R}(scale=ustrip(q_si), offset=0.0, basedim=dimension(q_si))
    return constructorof(Q)(one(T), dims)
end

#Conversions
for (type, _, _) in ABSTRACT_QUANTITY_TYPES
    @eval begin
        function Base.convert(::Type{Q}, q::UnionAbstractQuantity{<:Any,<:Dimensions}) where {T,Q<:$type{T,AffineDimensions}}
            return convert(with_type_parameters(Q, T, AffineDimensions{DEFAULT_DIM_BASE_TYPE}), q)
        end

        #Conversion of (AbstractQuantity){T,Dimensions{R}} to (AbstractQuantity){T,AffineDimensions{R}}
        function Base.convert(::Type{Q}, q::UnionAbstractQuantity{<:Any,<:Dimensions}) where {T,R,Q<:$type{T,AffineDimensions{R}}}
            dims = AffineDimensions{R}(scale=1, offset=0, basedim=dimension(q))
            return constructorof(Q)(convert(T, ustrip(q)), dims)
        end

        #Forced conversion of (AbstractQuantity){T,R<:AffineDimensions} to (AbstractQuantity){T,R<:Dimensions} (zero offset requirement overridden)
        function force_convert(::Type{Q}, q::UnionAbstractQuantity{<:Any,<:AbstractAffineDimensions}) where {T,D<:Dimensions,Q<:$type{T,D}}
            d = dimension(q)
            v = ustrip(q)*scale(d) + offset(d)
            return constructorof(Q)(convert(T, v), basedim(d))
        end

        #Conversion of (AbstractQuantity){T,R<:AffineDimensions} to (AbstractQuantity){T,R<:Dimensions}
        function Base.convert(::Type{Q}, q::UnionAbstractQuantity{<:Any,<:AbstractAffineDimensions}) where {T,D<:Dimensions,Q<:$type{T,D}}
            assert_no_offset(dimension(q))
            return force_convert(Q, q)
        end
    end
end

#Promotion rules
function Base.promote_rule(::Type{AffineDimensions{R1}}, ::Type{AffineDimensions{R2}}) where {R1,R2}
    return AffineDimensions{promote_type(R1,R2)}
end
function Base.promote_rule(::Type{AffineDimensions{R1}}, ::Type{Dimensions{R2}}) where {R1,R2}
    return Dimensions{promote_type(R1,R2)}
end
function Base.promote_rule(::Type{Dimensions{R1}}, ::Type{AffineDimensions{R2}}) where {R1,R2}
    return Dimensions{promote_type(R1,R2)}
end
function Base.promote_rule(::Type{SymbolicDimensions{R1}}, ::Type{AffineDimensions{R2}}) where {R1,R2}
    return Dimensions{promote_type(R1,R2)}
end
function Base.promote_rule(::Type{AffineDimensions{R1}}, ::Type{SymbolicDimensions{R2}}) where {R1,R2}
    return Dimensions{promote_type(R1,R2)}
end




# Conversions for Dimensions |> AffineDimenions =====================================================================================
"""
    uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractAffineDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions})

Convert a quantity `q` with base SI units to the affine units of `qout`, for `q` and `qout` with compatible units.
You can also use `|>` as a shorthand for `uconvert`
"""
function uconvert(qout::UnionAbstractQuantity{<:Any, <:AffineDimensions}, q::UnionAbstractQuantity{<:Any, <:Dimensions})
    @assert isone(ustrip(qout)) "You passed a quantity with a non-unit value to uconvert."
    dout = dimension(qout)
    dimension(q) == basedim(dout) || throw(DimensionError(q, qout))
    vout = (ustrip(q)-offset(dout))/scale(dout)
    return new_quantity(typeof(q), vout, dout)
end

function uconvert(qout::UnionAbstractQuantity{<:Any,<:AffineDimensions}, q::QuantityArray{<:Any,<:Any,<:Dimensions})
    @assert isone(ustrip(qout)) "You passed a quantity with a non-unit value to uconvert."
    dout = dimension(qout)
    dimension(q) == basedim(dout) || throw(DimensionError(q, qout))
    vout = (ustrip(q) .- offset(dout))./scale(dout)
    return QuantityArray(vout, dout, quantity_type(q))
end

# Conversions for AbstractAffineDimensions |> AbstractSymbolicDimensions =======================================================
function uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractSymbolicDimensions}, qin::AbstractQuantityOrArray{<:Any, <:AbstractAffineDimensions})
    uconvert(qout, siunits(qin))
end

# Conversions for AbstractSymbolicDimensions |> AbstractAffineDimensions   =======================================================
function uconvert(qout::UnionAbstractQuantity{<:Any,<:AbstractAffineDimensions}, qin::AbstractQuantityOrArray{<:Any, <:AbstractSymbolicDimensions})
    uconvert(qout, siunits(qin))
end

# Conversions for AbstractAffineDimensions |> AbstractAffineDimensions =======================================================
function uconvert(qout::UnionAbstractQuantity{<:Any, <:AbstractAffineDimensions}, qin::AbstractQuantityOrArray{<:Any, <:AbstractAffineDimensions})
    uconvert(qout, siunits(qin))
end

# Multiplication and division of AffineDimensions ===============================================================
function Base.:*(l::AffineDimensions, r::AffineDimensions) 
    assert_no_offset(l)
    assert_no_offset(r)
    return AffineDimensions(
        scale  = scale(l)*scale(r),
        offset = offset(l),
        basedim = basedim(l)*basedim(r)
    )
end

function Base.:/(l::AffineDimensions, r::AffineDimensions) 
    assert_no_offset(l)
    assert_no_offset(r)
    return AffineDimensions(
        scale  = scale(l)/scale(r),
        offset = offset(l),
        basedim = basedim(l)/basedim(r)
    )
end

# Exponentiation ===============================================================
function Base.:^(l::AffineDimensions{R}, r::Number) where {R}
    assert_no_offset(l)
    return AffineDimensions(
        scale   = scale(l)^r,
        offset  = offset(l),
        basedim = basedim(l)^tryrationalize(R, r)
    )
end

function Base.:^(l::AffineDimensions{R}, r::Integer) where {R}
    assert_no_offset(l)
    return AffineDimensions(
        scale   = scale(l)^r,
        offset  = offset(l),
        basedim = basedim(l)^tryrationalize(R, r)
    )
end

function Base.:inv(l::AffineDimensions{R}) where {R}
    assert_no_offset(l)
    return AffineDimensions(
        scale   = inv(scale(l)),
        offset  = offset(l),
        basedim = inv(basedim(l))
    )
end

# Operations on self-values ======================================================================================
function _no_offset_expand(q::Q) where {T, R, D<:AbstractAffineDimensions{R}, Q<:UnionAbstractQuantity{T,D}}
    return convert(with_type_parameters(Q, T, Dimensions{R}), q)
end

#Addition will return Quantity{T, Dimensions}
Base.:+(q1::UnionAbstractQuantity{<:Any,<:AffineDimensions}, q2::UnionAbstractQuantity{<:Any,<:AffineDimensions}) = _no_offset_expand(q1) + _no_offset_expand(q2)

#Subtraction will return Quantity{T, Dimensions}, in special cases, differences between offsetted AffineDimensions is allowed as offsets cancel out
function Base.:-(q1::UnionAbstractQuantity{<:Any,<:AffineDimensions}, q2::UnionAbstractQuantity{<:Any,<:AffineDimensions})
    if dimension(q1) == dimension(q2)
        return siunits(q1) - siunits(q2)
    else
        return _no_offset_expand(q1) - _no_offset_expand(q2)
    end
end

Base.:(==)(q1::UnionAbstractQuantity{<:Any, <:AffineDimensions}, q2::UnionAbstractQuantity{<:Any, <:AffineDimensions})   = (siunits(q1) == siunits(q2))
Base.:(==)(q1::UnionAbstractQuantity{<:Any, <:AffineDimensions}, q2::UnionAbstractQuantity{<:Any, <:AbstractDimensions}) = (siunits(q1) == siunits(q2))
Base.:(==)(q1::UnionAbstractQuantity{<:Any, <:AbstractDimensions}, q2::UnionAbstractQuantity{<:Any, <:AffineDimensions}) = (siunits(q1) == siunits(q2))
Base.:(≈)(q1::UnionAbstractQuantity{<:Any, <:AffineDimensions}, q2::UnionAbstractQuantity{<:Any, <:AffineDimensions})    = (siunits(q1) ≈ siunits(q2))
Base.:(≈)(q1::UnionAbstractQuantity{<:Any, <:AffineDimensions}, q2::UnionAbstractQuantity{<:Any, <:AbstractDimensions})  = (siunits(q1) ≈ siunits(q2))
Base.:(≈)(q1::UnionAbstractQuantity{<:Any, <:AbstractDimensions}, q2::UnionAbstractQuantity{<:Any, <:AffineDimensions})  = (siunits(q1) ≈ siunits(q2))
Base.:(==)(d1::AffineDimensions, d2::AffineDimensions) = (d1.scale==d2.scale) & (d1.offset==d2.offset) & (d1.basedim == d2.basedim)

# Units are stored using SymbolicDimensionsSingleton
const DEFAULT_AFFINE_QUANTITY_TYPE = with_type_parameters(DEFAULT_QUANTITY_TYPE, DEFAULT_VALUE_TYPE, AffineDimensions{DEFAULT_DIM_BASE_TYPE})

module AffineUnits

    using DispatchDoctor: @unstable

    import ..affine_unit
    import ..scale
    import ..offset
    import ..basedim
    import ..dimension
    import ..ustrip
    import ..constructorof
    import ..DEFAULT_AFFINE_QUANTITY_TYPE
    import ..DEFAULT_DIM_TYPE
    import ..DEFAULT_VALUE_TYPE
    import ..Units: UNIT_SYMBOLS, UNIT_VALUES
    import ..Constants: CONSTANT_SYMBOLS, CONSTANT_VALUES
    import ..Constants
    import ..Quantity
    import ..INDEX_TYPE
    import ..AbstractDimensions
    import ..AffineDimensions
    import ..UnionAbstractQuantity

    import ..DEFAULT_DIM_BASE_TYPE
    import ..WriteOnceReadMany
    import ..SymbolicUnits.as_quantity

    #Constants are not imported
    const AFFINE_UNIT_SYMBOLS = WriteOnceReadMany([UNIT_SYMBOLS...])
    const AFFINE_UNIT_VALUES  = WriteOnceReadMany(affine_unit.([UNIT_VALUES...]))
    const AFFINE_UNIT_MAPPING = WriteOnceReadMany(Dict(s => INDEX_TYPE(i) for (i, s) in enumerate(AFFINE_UNIT_SYMBOLS)))

    # Used for registering units in current module
    function update_external_affine_unit(name::Symbol, q::UnionAbstractQuantity{<:Any,<:AffineDimensions{R}}) where R
        ind = get(AFFINE_UNIT_MAPPING, name, INDEX_TYPE(0))
        if !iszero(ind)
            @warn "unit $(name) already exists, skipping"
            return nothing
        end

        #Extract original dimensions
        dims  = dimension(q)

        #Add "name" to the symbol to make it display
        d_sym = AffineDimensions{DEFAULT_DIM_BASE_TYPE}(
            scale = scale(dims),
            offset = offset(dims),
            basedim = basedim(dims),
            symbol = name
        )

        #Reconstruct the quantity with the new name
        q_sym = constructorof(DEFAULT_AFFINE_QUANTITY_TYPE)(ustrip(q), d_sym)

        push!(AFFINE_UNIT_SYMBOLS, name)
        push!(AFFINE_UNIT_VALUES, q_sym)
        AFFINE_UNIT_MAPPING[name] = lastindex(AFFINE_UNIT_SYMBOLS)
        return nothing
    end
    
    update_external_affine_unit(name::Symbol, q::UnionAbstractQuantity) = update_external_affine_unit(name, affine_unit(q))
    update_external_affine_unit(name::Symbol, d::AbstractDimensions)    = update_external_affine_unit(name, Quantity(DEFAULT_VALUE_TYPE(1.0), d))
    function update_external_affine_unit(d::AffineDimensions)  
        d.symbol != :nothing || error("Cannot register affine dimension if symbol is :nothing")
        return update_external_affine_unit(d.symbol, d)
    end

    """
        aff_uparse(s::AbstractString)

    Parse a string containing an expression of units and return the
    corresponding `Quantity` object with `Float64` value. 
    However, unlike the regular `u"..."` macro, this macro uses
    `AffineDimensions` for the dimension type, which can represent a greater
    number of units, but much more limited functionality with calculations. 
    For example, `aff_uparse("km/s^2")` would be parsed to
    `Quantity(1.0, AffineDimensions(scale=1000.0, offset=0.0, basedim=Dimensions(length=1, time=-2)))`.
    """
    function aff_uparse(s::AbstractString)
        ex = map_to_scope(Meta.parse(s))
        ex = :($as_quantity($ex))
        return eval(ex)::DEFAULT_AFFINE_QUANTITY_TYPE
    end

    as_quantity(q::DEFAULT_AFFINE_QUANTITY_TYPE) = q

    """
        ua"[unit expression]"

    Parse a string containing an expression of units and return the
    corresponding `Quantity` object with `Float64` value. 
    However, unlike the regular `u"..."` macro, this macro uses
    `AffineDimensions` for the dimension type, which can represent a greater
    number of units, but much more limited functionality with calculations. 
    For example, `ua"km/s^2"` would be parsed to
    `Quantity(1.0, AffineDimensions(scale=1000.0, offset=0.0, basedim=Dimensions(length=1, time=-2)))`.
    """
    macro ua_str(s)
        ex = map_to_scope(Meta.parse(s))
        ex = :($as_quantity($ex))
        return esc(ex)
    end

    @unstable function map_to_scope(ex::Expr)
        if !(ex.head == :call)
            throw(ArgumentError("Unexpected expression: $ex. Only `:call` is expected."))
        end
        if ex.head == :call
            ex.args[2:end] = map(map_to_scope, ex.args[2:end])
            return ex
        end
    end

    function map_to_scope(sym::Symbol)
        if sym in AFFINE_UNIT_SYMBOLS
            return lookup_unit(sym)
        else
            throw(ArgumentError("Symbol $sym not found in `AffineUnits`."))
        end
    end

    function map_to_scope(ex)
        return ex
    end

    function lookup_unit(ex::Symbol)
        i = findfirst(==(ex), AFFINE_UNIT_SYMBOLS)::Int
        return AFFINE_UNIT_VALUES[i]
    end

    #Register Celsius and Fahrenheit (the most commonly used affine units)
    begin
        K  = Quantity(1.0, temperature=1)
        °C = Quantity(1.0, AffineDimensions(scale=1.0, offset=273.15*K, basedim=K, symbol=:°C))
        °F = Quantity(1.0, AffineDimensions(scale=5/9, offset=(-160/9)°C, basedim=°C, symbol=:°F))
        update_external_affine_unit(dimension(°C))
        update_external_affine_unit(:degC, dimension(°C))
        update_external_affine_unit(dimension(°F))
        update_external_affine_unit(:degF, dimension(°F))
    end
    
end



import .AffineUnits: aff_uparse, update_external_affine_unit

"""
    ua"[unit expression]"

    Parse a string containing an expression of units and return the
    corresponding `Quantity` object with `Float64` value. 
    However, unlike the regular `u"..."` macro, this macro uses
    `AffineDimensions` for the dimension type, which can represent a greater
    number of units, but supports a much smaller set of operations. It is
    adviced to convert AffineDimensions to regular are symbolic dimensions
    as soon as possible. 
    For example, `ua"km/s^2"` would be parsed to
    `Quantity(1.0, AffineDimensions(scale=1000.0, offset=0.0, basedim=Dimensions(length=1, time=-2)))`.
"""
macro ua_str(s)
    ex = AffineUnits.map_to_scope(Meta.parse(s))
    ex = :($AffineUnits.as_quantity($ex))
    return esc(ex)
end







