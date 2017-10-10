function differentiate(F::Vector{FP.Polynomial{T}}) where {T<:Number}
    [FP.differentiate(f, i) for f in F, i=1:FP.nvariables.(F[1])]
end
