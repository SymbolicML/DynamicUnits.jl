module DynamicQuantities

export Units, Constants, SymbolicUnits, SymbolicConstants
export AbstractQuantity, AbstractGenericQuantity, AbstractRealQuantity, UnionAbstractQuantity
export Quantity, GenericQuantity, RealQuantity
export FixedRational
export AbstractDimensions, Dimensions, NoDims
export AbstractSymbolicDimensions, SymbolicDimensions, SymbolicDimensionsSingleton
export AbstractAffineDimensions, AffineDimensions
export QuantityArray
export DimensionError
export ustrip, dimension
export ulength, umass, utime, ucurrent, utemperature, uluminosity, uamount
export uparse, @u_str, sym_uparse, @us_str, aff_uparse, @ua_str, uexpand, uconvert, @register_unit, @register_affine_unit

# Deprecated:
export expand_units

using DispatchDoctor: @stable

@stable default_mode="disable" begin
    include("internal_utils.jl")
    include("fixed_rational.jl")
    include("write_once_read_many.jl")
    include("types.jl")
    include("utils.jl")
    include("math.jl")
    include("arrays.jl")
    include("units.jl")
    include("constants.jl")
    include("uparse.jl")
    include("symbolic_dimensions.jl")
    include("complex.jl")
    include("register_units.jl")
    include("disambiguities.jl")
    include("affine_dimensions.jl")

    include("deprecated.jl")
end


import .Units
import .Constants
import .UnitsParse: uparse, @u_str

using .Units: UNIT_SYMBOLS

# Copy all units to top level:
let _units_import_expr = :(using .Units: m, g)
    append!(
        _units_import_expr.args[1].args,
        Expr(:(.), s) for s in UNIT_SYMBOLS if s ∉ (:m, :g)
    )
    eval(_units_import_expr)
end

end
