#Simulates svd for bigfloat matrices by converting to Float64/Complex128.
Base.LinAlg.svdfact!(A::Matrix{BigFloat}; thin=true) = svdfact!(convert.(Float64, A); thin=true)
Base.LinAlg.svdfact!(A::Array{Complex{BigFloat},2}; thin=true) = svdfact!(convert.(Complex128, A); thin=true)
