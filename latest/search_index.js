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
    "text": "Homotopy.jl is a package for constructing (polynomial) homotopies H(xt).Each implemented homotopy has the same Interface so that you can switch easily between different homotopy types. Based on this interface there are also some convenient higher level constructs provided, e.g. the construction of a total degree system and its start solutions."
},

{
    "location": "index.html#Example-1",
    "page": "Introduction",
    "title": "Example",
    "category": "section",
    "text": "using Homotopy\n# we use an MultivariatePolynomials implementation to construct the homotopy.\nimport DynamicPolynomials: @polyvar\n\n@polyvar x y z\n\nH = StraightLineHomotopy([x + y^3, x^2*y-2y], [x^3+2, y^3+2])\n# H is now StraightLineHomotopy{Int64},\n# but let's assume our algorithm uses Complex128, to avoid unnecessary conversions\n# it would be better to make\nH = StraightLineHomotopy{Complex128}([x + y^3, x^2*y-2y], [x^3+2, y^3+2])\n\n# we can now evaluate H\nevaluate(H, rand(Complex128, 2), 0.42)\n# or alternatively\nH(rand(Complex128, 2), 0.42)"
},

{
    "location": "index.html#Homotopy.StraightLineHomotopy",
    "page": "Introduction",
    "title": "Homotopy.StraightLineHomotopy",
    "category": "Type",
    "text": "StraightLineHomotopy(start, target)\n\nConstruct the homotopy t * start + (1-t) * target.\n\nstart and target have to match and to be one of the following\n\nVector{<:MP.AbstractPolynomial} where MP is MultivariatePolynomials\nMP.AbstractPolynomial\nVector{<:FP.Polynomial} where FP is FixedPolynomials\n\nStraightLineHomotopy{T}(start, target)\n\nYou can also force a specific coefficient type T.\n\n\n\n"
},

{
    "location": "index.html#Homotopy.GammaTrickHomotopy",
    "page": "Introduction",
    "title": "Homotopy.GammaTrickHomotopy",
    "category": "Type",
    "text": "GammaTrickHomotopy(start, target, [γ])\n\nConstruct the homotopy (1-t)⋅target + t⋅γ⋅start. If γ is not supplied it will be drawn randomly and uniformly from the complex unit circle.\n\nstart and target have to match and to be one of the following\n\nVector{<:MP.AbstractPolynomial} where MP is MultivariatePolynomials\nMP.AbstractPolynomial\nVector{<:FP.Polynomial} where FP is FixedPolynomials\n\nGammaTrickHomotopy(start, target, seed::Int)\n\nYou can also supply a seed for the RNG which is used to draw γ, i.e. subsequent invocations with the same seed will produce the same γ.\n\nGammaTrickHomotopy{T}(start, target, [γ])\nGammaTrickHomotopy{T}(start, target, seed::Int)\n\nYou can also force a specific coefficient type T.\n\n\n\n"
},

{
    "location": "index.html#Homotopies-1",
    "page": "Introduction",
    "title": "Homotopies",
    "category": "section",
    "text": "The following homotopies are implementedStraightLineHomotopy\nGammaTrickHomotopy"
},

{
    "location": "interface.html#",
    "page": "Interface",
    "title": "Interface",
    "category": "page",
    "text": ""
},

{
    "location": "interface.html#Interface-1",
    "page": "Interface",
    "title": "Interface",
    "category": "section",
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
    "location": "interface.html#Evaluation-1",
    "page": "Interface",
    "title": "Evaluation",
    "category": "section",
    "text": "evaluate\nevaluate!"
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
    "location": "interface.html#Differentiation-1",
    "page": "Interface",
    "title": "Differentiation",
    "category": "section",
    "text": "jacobian\njacobian!\ndt\ndt!"
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
    "location": "interface.html#Homogenization-1",
    "page": "Interface",
    "title": "Homogenization",
    "category": "section",
    "text": "homogenize\ndehomogenize\nishomogenized\nishomogenous"
},

{
    "location": "interface.html#Homotopy.nvariables",
    "page": "Interface",
    "title": "Homotopy.nvariables",
    "category": "Function",
    "text": "nvariables(H::AbstractHomotopy)\n\nThe number of variables which H expects as input, i.e. to evaluate H(x,t) x has to be a vector of length nvariables(H).\n\n\n\n"
},

{
    "location": "interface.html#Misc-1",
    "page": "Interface",
    "title": "Misc",
    "category": "section",
    "text": "nvariables"
},

{
    "location": "higherlevelconstructs.html#",
    "page": "Higher level constructs",
    "title": "Higher level constructs",
    "category": "page",
    "text": ""
},

{
    "location": "higherlevelconstructs.html#higherlevelconstructs-1",
    "page": "Higher level constructs",
    "title": "Higher level constructs",
    "category": "section",
    "text": ""
},

{
    "location": "higherlevelconstructs.html#Homotopy.totaldegree",
    "page": "Higher level constructs",
    "title": "Homotopy.totaldegree",
    "category": "Function",
    "text": "totaldegree(H::Type{AbstractHomotopy}, F, [unitroots=false])\n\nConstruct a  total degree homotopy of type H with F and an iterator of its solutions. This is the homotopy with start system\n\nbeginalign*\n    z_1^d_1 - b_1\n    z_1^d_2 - b_2\n    vdots \n    z_n^d_n - b_n\nendalign*\n\nand target system F, where d_i is the degree of the i-th polynomial of F. If unitroots=true then b_i=1 otherwise b_i is a random complex number (with real and imaginary part in the unit interval).\n\nExample\n\nH, startsolutions = totaldegree(StraightLineHomotopy{Complex128}, [x^2+y+1, x^3*y-2])\n\n\n\n"
},

{
    "location": "higherlevelconstructs.html#Homotopy.TotalDegreeSolutionIterator",
    "page": "Higher level constructs",
    "title": "Homotopy.TotalDegreeSolutionIterator",
    "category": "Type",
    "text": "TotalDegreeSolutionIterator(degrees, b)\n\nGiven the Vectors degrees and b TotalDegreeSolutionIterator enumerates all solutions of the system\n\nbeginalign*\n    z_1^d_1 - b_1 = 0 \n    z_1^d_2 - b_2 = 0 \n    vdots \n    z_n^d_n - b_n = 0 \nendalign*\n\nwhere d_i is degrees[i] and b_i is b[i].\n\n\n\n"
},

{
    "location": "higherlevelconstructs.html#Homotopy.totaldegree_startsystem",
    "page": "Higher level constructs",
    "title": "Homotopy.totaldegree_startsystem",
    "category": "Function",
    "text": "totaldegree_startsystem(F::Vector{FP.Polynomial{<:Complex}}, [unit_roots=false])\n\nReturn the system\n\nbeginalign*\n    z_1^d_1 - b_1\n    z_1^d_2 - b_2\n    vdots \n    z_n^d_n - b_n\nendalign*\n\nwhere d_i is the degree of the i-th polynomial of F and an iterator of its solutions. If unitroots=true then b_i=1 otherwise b_i is a random complex number (with real and imaginary part in the unit interval).\n\n\n\n"
},

{
    "location": "higherlevelconstructs.html#Total-degree-homotopy-1",
    "page": "Higher level constructs",
    "title": "Total degree homotopy",
    "category": "section",
    "text": "totaldegree\nTotalDegreeSolutionIterator\ntotaldegree_startsystem"
},

{
    "location": "higherlevelconstructs.html#Homotopy.randomhomotopy",
    "page": "Higher level constructs",
    "title": "Homotopy.randomhomotopy",
    "category": "Function",
    "text": "randomhomotopy(::Type{AbstractHomotopy{T}}, size::Int; kwargs...)\n\nCreate a total degree homotopy where the target system is a randomsystem(T, size, size; kwargs...).\n\nExample\n\njulia> H, solutions = randomhomotopy(StraightLineHomotopy{Complex128}, 2, mindegree=3, maxdegree=6);\njulia> length(H)\n3\njulia> nvariables(H)\n3\n\n\n\n"
},

{
    "location": "higherlevelconstructs.html#Homotopy.randomsystem",
    "page": "Higher level constructs",
    "title": "Homotopy.randomsystem",
    "category": "Function",
    "text": "randomsystem([T=Complex128,] nequations::Int, nvars::Int; mindegree=0, maxdegree=5, rng=Base.Random.GLOBAL_RNG)\n\nCreates a random polynomial system of nequations equations with nvars variables (named x_1, ...x_nvars). Each polynomial has a total degree uniformly drawn from mindegree maxdegree. The coefficients are drawn independently from the given rng.\n\nrandomsystem([T=Complex128,] degrees::Vector{Int}, variables::Vector{Symbol}; rng=N(0,1))\n\nCreate a random polynomial system with the given degrees and variables.\n\n\n\n"
},

{
    "location": "higherlevelconstructs.html#Random-homotopies-1",
    "page": "Higher level constructs",
    "title": "Random homotopies",
    "category": "section",
    "text": "randomhomotopy\nrandomsystem"
},

]}
