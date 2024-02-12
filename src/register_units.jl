import .Units: UNIT_MAPPING, UNIT_SYMBOLS, UNIT_VALUES, _lazy_register_unit
import .SymbolicUnits: update_external_symbolic_unit_value

# Update the unit collections
const UNIT_UPDATE_LOCK = Threads.SpinLock()

function update_all_values(name_symbol, unit)
    lock(UNIT_UPDATE_LOCK) do
        push!(ALL_SYMBOLS, name_symbol)
        push!(ALL_VALUES, unit)
        i = lastindex(ALL_VALUES)
        ALL_MAPPING[name_symbol] = i
        UNIT_MAPPING[name_symbol] = i
        update_external_symbolic_unit_value(name_symbol)
    end
end

"""
    @register_unit symbol value

Register a new unit under the given symbol to have
a particular value.
"""
macro register_unit(symbol, value)
    return esc(_register_unit(symbol, value))
end

function _register_unit(name::Symbol, value)
    name_symbol = Meta.quot(name)
    index = get(ALL_MAPPING, name, INDEX_TYPE(0))
    if iszero(index)
        reg_expr = _lazy_register_unit(name, value)
        push!(reg_expr.args, quote
            $update_all_values($name_symbol, $value)
            nothing
        end)
        return reg_expr
    else
        unit = ALL_VALUES[index]
        # When a utility function to expand `value` to its final form becomes
        # available, enable the following check. This will avoid throwing an error
        # if user is trying to register an existing unit with matching values.
        # unit.value != value && throw("Unit $name is already defined as $unit")
        throw("Unit `$name` is already defined as `$unit`")
    end
end
