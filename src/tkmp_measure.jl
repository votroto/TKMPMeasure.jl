using MultivariatePolynomials
using MultivariateMoments
using SemialgebraicSets
using LinearAlgebra
using JuMP

struct TKMPMeasureBuilder
	maxdegree::Int
	on
	vars
end

mutable struct TKMPMeasure
	var_map
	vars
	maxdegree::Int
end

struct NonnegativeOnKConstraintBuilder
	m::TKMPMeasure
	k::AbstractSemialgebraicSet
end

# -----------------------------------------------------------------------------
# utils
# -----------------------------------------------------------------------------

SemialgebraicSets.inequalities(::AbstractSemialgebraicSet) = []

function coeffs(p::AbstractPolynomialLike, vars)
	deg = maxdegree(p)
	ms = monomials(vars, 0:deg) |> reverse
	cs = MultivariatePolynomials.coefficient.(Ref(p), ms, Ref(vars))

	cs, ms
end

function integrate(p::AbstractPolynomialLike, d::TKMPMeasure)
	trimvec(v::AbstractVector, t) = v[begin:min(t, length(v))]

	l = length(d.var_map.data)
	cs, ms = trimvec.(coeffs(p, d.vars), l)
	us = map(m -> d.var_map[m], ms)

	dot(cs, us)
end

integrate(p::Number, d::TKMPMeasure) = integrate(p * only(monomials(d.vars, 0)), d)
integrate(p::AbstractPolynomialLike, d::AbstractVector{TKMPMeasure}) = reduce(integrate, d; init=p)
∫ = integrate

function _moment_matrix(m::TKMPMeasure)
	xd = monomials(m.vars, 0:m.maxdegree ÷ 2) |> reverse
	Md = xd * xd'
	integrate.(Md, Ref(m))
end

function localizing_matrix(m::TKMPMeasure, g::AbstractPolynomialLike)
	round_up_even(x) = cld(x,2)*2

	tg = round_up_even(maxdegree(g))
	xd = monomials(m.vars, 0:(m.maxdegree ÷ 2 - tg ÷ 2)) |> reverse
	Ld = xd * xd' * g
	integrate.(Ld, Ref(m))
end

function JuMP.value(m::TKMPMeasure)
	polyeval(p::AbstractPolynomialLike) = value(p)((variables(p) .=> 1)...)

	b = monomials(m.vars, 0:m.maxdegree ÷ 2) |> reverse
	m = polyeval.(_moment_matrix(m))

	extractatoms(MomentMatrix(m, b), 1e-3)
end

# -----------------------------------------------------------------------------
# JuMP Variables
# -----------------------------------------------------------------------------

function JuMP.build_variable(_err::Function, info::VariableInfo, ::Type{TKMPMeasure}; maxdegree, on)
	TKMPMeasureBuilder(maxdegree, on, variables(on))
end

function JuMP.add_variable(model::JuMP.Model, b::TKMPMeasureBuilder, name::String)
	ms = monomials(b.vars, 0:b.maxdegree) |> reverse
	var_map = @variable(model, [ms]; base_name=name)

	μ = TKMPMeasure(var_map, b.vars, b.maxdegree)

	@constraint(model, _moment_matrix(μ) in PSDCone(), base_name=string(name," moment matrix"))
	@constraint(model, [i in inequalities(b.on)], localizing_matrix(μ, i) in PSDCone(), base_name=string(name," localizing matrix inequalities"))
	@constraint(model, [i in equalities(b.on)], localizing_matrix(μ, i) .== 0, base_name=string(name," localizing matrix equalities"))

	μ
end
