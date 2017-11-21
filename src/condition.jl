export κ, κ_norm, μ_norm, kappa, kappa_norm, mu_norm


"""
    κ(H, z, t, cfg)

Computes the condition number of H at (z, t) (with config cfg). See Condition^[1] for details

Proposition 16.10: κ(f,z) := ‖f‖ ‖ Df(z)^† diag(‖ z ‖^{d_i-1}) ‖

[1]: Condition, Bürgisser and Cucker
"""


function κ(H::AbstractHomotopy, z::Vector, t::Float64, cfg=PolynomialHomotopyConfig(H))
    weyl_norm=weylnorm(H)
    f = weyl_norm(t)
    D = jacobian(H,z,t,cfg)
    norm_z = norm(z)

    M = diagm(map(d_i -> norm_z ^ (1-d_i), FP.degree.(H.start)))

    try
        _, s, _ = svd(M*D)
        σ = s[end]

        return real(f * inv(σ))
    catch
        Inf
    end
end



kappa(H::AbstractHomotopy, z::Vector, t::Float64, cfg=PolynomialHomotopyConfig(H))=κ(H, z, t, cfg)



"""
    κ_norm(H, z, t, cfg)

Computes the condition number of H at (z, t) (with config cfg). See Condition^[1] for details

Eq. (16.11): κ_norm(f,z) := ‖f‖ ‖ Df(z)^† diag(√{d_i}‖ z ‖^{d_i-1}) ‖

[1]: Condition, Bürgisser and Cucker
"""


function κ_norm(H::AbstractHomotopy, z::Vector, t::Float64, cfg=PolynomialHomotopyConfig(H))
    weyl_norm=weylnorm(H)
    f = weyl_norm(t)
    D = jacobian(H,z,t,cfg)
    norm_z = norm(z)

    M = diagm(map(d_i -> inv(sqrt(d_i)) * norm_z ^ (1-d_i), FP.degree.(H.start)))

    _, s, _ = svd(M*D)
    σ = s[end]

    real(f * inv(σ))
end

kappa_norm(H::AbstractHomotopy, z::Vector, t::Float64, cfg=PolynomialHomotopyConfig(H))=κ_norm(H, z, t, cfg)


"""
    μ_norm(H, z, t, cfg)

Computes the condition number of H at (z, t) (with config cfg). See Condition^[1] for details

Definition 16.43: μ_norm(f,z) := ‖f‖ ‖ (Df(z)-(T_z))^{-1} diag(√{d_i}‖ z ‖^{d_i-1}) ‖

[1]: Condition, Bürgisser and Cucker
"""


function μ_norm(H::AbstractHomotopy, z::Vector, t::Float64, cfg=PolynomialHomotopyConfig(H))
    if ishomogenized(H) || ishomogenous(H)

        weyl_norm=weylnorm(H)
        f = weyl_norm(t)
        D = jacobian(H,z,t,cfg)
        n=length(z)
        norm_z = norm(z)
        P = qrfact(reshape(z,n,1))
        Q = P[:Q][:,2:end]

        M = diagm(map(d_i -> inv(sqrt(d_i)) * norm_z ^ (1-d_i), FP.degree.(H.start)))

        _, s, _ = svd(M*D*Q)
        σ = last(s)

        real(f * inv(σ))
    else
        κ_norm(H, z, t, cfg)
    end
end

mu_norm(H::AbstractHomotopy, z::Vector, t::Float64, cfg=PolynomialHomotopyConfig(H))=μ_norm(H, z, t, cfg)
