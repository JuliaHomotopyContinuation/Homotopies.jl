var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "Homotopy.jl is a package for constructing (polynomial) homotopies H(xt).Each implemented homotopy has the same Interface so that you can switch easily between different homotopy types."
},

{
    "location": "index.html#Homotopy.StraightLineHomotopy",
    "page": "Introduction",
    "title": "Homotopy.StraightLineHomotopy",
    "category": "Type",
    "text": "StraightLineHomotopy(start, target)\n\nConstruct the homotopy t * start + (1-t) * target.\n\nstart and target have to match and to be one of the following\n\nVector{<:MP.AbstractPolynomial} where MP is MultivariatePolynomials\nMP.AbstractPolynomial\nVector{<:FP.Polynomial} where FP is FixedPolynomials\nStraightLineHomotopy{T}(start, target)\n\nYou can also force a specific coefficient type T.\n\n\n\n"
},

{
    "location": "index.html#Homotopy.GammaTrickHomotopy",
    "page": "Introduction",
    "title": "Homotopy.GammaTrickHomotopy",
    "category": "Type",
    "text": "GammaTrickHomotopy(start, target, [γ])\n\nConstruct the homotopy (1-t)⋅target + t⋅γ⋅start. If γ is not supplied it will be drawn randomly and uniformly from the complex unit circle.\n\nstart and target have to match and to be one of the following\n\nVector{<:MP.AbstractPolynomial} where MP is MultivariatePolynomials\nMP.AbstractPolynomial\nVector{<:FP.Polynomial} where FP is FixedPolynomials\nGammaTrickHomotopy(start, target, seed::Int)\n\nYou can also supply a seed for the RNG which is used to draw γ, i.e. subsequent invocations with the same seed will produce the same γ.\n\nGammaTrickHomotopy{T}(start, target, [γ])\nGammaTrickHomotopy{T}(start, target, seed::Int)\n\nYou can also force a specific coefficient type T.\n\n\n\n"
},

{
    "location": "index.html#Homotopies-1",
    "page": "Introduction",
    "title": "Homotopies",
    "category": "section",
    "text": "The following homotopies are implementedStraightLineHomotopy\nGammaTrickHomotopy"
},

{
    "location": "index.html#Example-1",
    "page": "Introduction",
    "title": "Example",
    "category": "section",
    "text": "using Homotopy\n# we use an MultivariatePolynomials implementation to construct the homotopy.\nimport DynamicPolynomials: @polyvar\n\n@polyvar x y z\n\nH = StraightLineHomotopy([x + y^3, x^2*y-2y], [x^3+2, y^3+2])\n# H is now StraightLineHomotopy{Int64},\n# but let's assume our algorithm uses Complex128, to avoid unnecessary conversions\n# it would be better to make\nH = StraightLineHomotopy{Complex128}([x + y^3, x^2*y-2y], [x^3+2, y^3+2])\n\n# we can now evaluate H\nevaluate(H, rand(Complex128, 2), 0.42)\n# or alternatively\nH(rand(Complex128, 2), 0.42)"
},

{
    "location": "interface.html#",
    "page": "Interface",
    "title": "Interface",
    "category": "page",
    "text": ""
},

{
    "location": "interface.html#Homotopy.evaluate",
    "page": "Interface",
    "title": "Homotopy.evaluate",
    "category": "Function",
    "text": "evaluate(H::AbstractHomotopy, x, t)\n\nEvaluate the homotopy H at x to time t, i.e. H(xt)\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.evaluate!",
    "page": "Interface",
    "title": "Homotopy.evaluate!",
    "category": "Function",
    "text": "evaluate!(u::Vector, H::AbstractHomotopy, x, t)\n\nEvaluate the homotopy H at x to time t, i.e. H(xt), and store the result in u. Use this instead of evaluate to avoid allocations.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.jacobian",
    "page": "Interface",
    "title": "Homotopy.jacobian",
    "category": "Function",
    "text": "jacobian(H::AbstractHomotopy)\n\nCompute an evaluation function (x, t) -> J_H(x,t) of the jacobian J_H of the homotopy H. The jacobian is constructed w.r.t. x, i.e. it doesn't contain the partial derivatives w.r.t. t.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.jacobian!",
    "page": "Interface",
    "title": "Homotopy.jacobian!",
    "category": "Function",
    "text": "jacobian!(H::AbstractHomotopy)\n\nCompute an inplace evaluation function (u, x, t) -> u := J_H(x,t) of the jacobian J_H of the homotopy H. Use this instead of jacobian to avoid allocations.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.dt",
    "page": "Interface",
    "title": "Homotopy.dt",
    "category": "Function",
    "text": "dt(H::AbstractHomotopy)\n\nCompute an evaluation function (x, t) -> ∂H∂t(x,t) of the partial derivative fracHt of the homotopy H.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.dt!",
    "page": "Interface",
    "title": "Homotopy.dt!",
    "category": "Function",
    "text": "dt!(H::AbstractHomotopy)\n\nCompute an inplace evaluation function (u, x, t) -> u := ∂H∂t(x,t) of the partial derivative fracHt of the homotopy H. Use this instead of dt to avoid allocations.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.nvariables",
    "page": "Interface",
    "title": "Homotopy.nvariables",
    "category": "Function",
    "text": "nvariables(H::AbstractHomotopy)\n\nThe number of variables which H expects as input, i.e. to evaluate H(x,t) x has to be a vector of length nvariables(H).\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.homogenize",
    "page": "Interface",
    "title": "Homotopy.homogenize",
    "category": "Function",
    "text": "homogenize(H::AbstractHomotopy)\n\nHomogenize the homotopy H. This adds an additional variable. If H is already homogenized, this is the identity.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.dehomogenize",
    "page": "Interface",
    "title": "Homotopy.dehomogenize",
    "category": "Function",
    "text": "dehomogenize(H::AbstractHomotopy)\n\nDehomogenize the homotopy H. This removes the first variable. If H is not homogenized, this is the identity.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.ishomogenized",
    "page": "Interface",
    "title": "Homotopy.ishomogenized",
    "category": "Function",
    "text": "ishomogenized(H::AbstractHomotopy)\n\nCheck whether the homotopy H was homogenized.\n\n\n\n"
},

{
    "location": "interface.html#Homotopy.ishomogenous",
    "page": "Interface",
    "title": "Homotopy.ishomogenous",
    "category": "Function",
    "text": "ishomogenous(H::AbstractHomotopy)\n\nCheck whether the homotopy H is homogenous. This does not imply that H was homogenized.\n\n\n\n"
},

{
    "location": "interface.html#Interface-1",
    "page": "Interface",
    "title": "Interface",
    "category": "section",
    "text": "evaluate\nevaluate!\njacobian\njacobian!\ndt\ndt!\nnvariables\nhomogenize\ndehomogenize\nishomogenized\nishomogenous"
},

]}
