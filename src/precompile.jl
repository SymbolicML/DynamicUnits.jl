import PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    types = (Float32, Float64)
    same_dim_bin_ops = (+, -, ==, <=, >=, isapprox, min, max)
    diff_dim_bin_ops = (*, /)
    una_ops = (
        sqrt, inv, cbrt, abs, abs2, -, x -> x^3, x -> x^5.2, iszero,
        one, isfinite, oneunit, zero, x -> show(devnull, x)
    )
    @compile_workload begin
        for T in types
            x = Quantity{T}(0.5, length=1, time=-1)
            y = Quantity(T(0.2), mass=1, time=-1, current=-2, temperature=3, luminosity=4, amount=5)
            z = Quantity(T(1.5), SymbolicDimensions(m=1, kg=-2, fs=-2))
            x_qa = QuantityArray([x, 2*x, 3*x])

            uexpand(z)
            uconvert(Quantity(1.0, SymbolicDimensions(cm=1, ps=-1)))(x)

            for op in same_dim_bin_ops
                op(x, 2 * x)
                op(z, 2 * z)
            end
            for op in diff_dim_bin_ops
                op(x, 2 * y)
                op(z, 2 * z)
                
                op(x, T(0.5))
                op(T(0.5), x)
                op(z, T(0.5))
                op(T(0.5), z)
            end
            for op in una_ops
                op(x)
                op(y)
                op(z)
                op.(x_qa)
            end
            show(devnull, x_qa)
        end
    end
end
