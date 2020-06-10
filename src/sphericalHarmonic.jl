#normalization factor
abstract type YLMNorm end

struct Schmidt{T} <: YLMNorm
end

struct Laplace{T} <: YLMNorm
end

struct Ylm{T}; end
struct Rlm{T}; end
struct Rlylm{T}; end

function ylmKCoefficient(N::Schmidt{T},l::Int64, m::Int64) where T

	k = one(T)
	for i in (l-m+1):(l+m)
	  k *= i
	end
	return sqrt(1/k)
end

function ylmKCoefficient(N::Laplace{T}, l::Int64, m::Int64) where T

  k = one(T)
  for i in (l-m+1):(l+m)
    k *= i
  end

  return sqrt((2*l+1) / (4*T(pi)*k))
end


function ylmCosSinPolynomial(m::Int64, x::Variable, y::Variable)

  sum = zero(x*y)
  for j in 0:div(m,2)
    sum += ((-1)^j)*binomial(m, 2j)*(y^(2j))*(x^(m-2j))
  end
  return sum
end

function ylmSinSinPolynomial(m::Int64, x::Variable, y::Variable)

  sum = zero(x*y)
  for j in 0:div((m-1),2)
    sum += ((-1)^j)*binomial(m, 2j + 1)*(y^(2j + 1))*(x^(m-2j-1))
  end
  return sum
end

"""
    ylm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable)
*Description:*  Calculation of the spherical harmonic for a given order (l,m) in Cartesian coordinates\\

*Input:*  `l`       - Degree of the spherical harmonic\\
          `m`       - Order of the spherical harmonic\\
          `x, y, z` - Cartesian coordinates\\

*Output:*  Spherical harmonic polynomial
"""
function Ylm{T}(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable; norm::YLMNorm=Laplace{T}()) where T

  if abs(m) > l
    throw(DomainError(m,"-l <= m <= l expected, but m = $m and l = $l."))
  end

  p = (z^2 - 1)^l

  for i = 1:l+abs(m)
    c = i <= l ? 1/(2one(T)*i) : one(T)
    p = c*differentiate(p, z)
  end

  if m > 0
    return sqrt(2one(T))*ylmKCoefficient(norm, l, m)*ylmCosSinPolynomial(m,x,y)*p
  elseif m < 0
    return sqrt(2one(T))*ylmKCoefficient(norm, l, abs(m))*ylmSinSinPolynomial(abs(m),x,y)*p
  else
    return ylmKCoefficient(norm, l, 0)*p
  end
end

# multiplying r^l*ylm(x,y,z)
function Rlylm{T}(l::Int64, m::Int, x::Variable, y::Variable, z::Variable; kwargs...) where T
	p = Ylm{T}(l,m,x,y,z;kwargs...)
	tout = []
	# Zerlegung des Polynoms in Terme:
	for t in terms(p)
		deg = degree(monomial(t)) # Gibt den gesamten Grad des Monoms an
		degR = l-deg # durch das Kürzen ergibt sich ein Grad von l-deg fuer r
		push!(tout,(x^2+y^2+z^2)^Int(degR/2)*t) # r² wird durch x²+y²+z² ersetzt
	end

	return polynomial(tout)
end

# solid harmonics
function Rlm{T}(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable; norm::YLMNorm=Laplace{T}()) where T
	rlm = Rlylm{T}(l,m,x,y,z; norm=norm)
	rlm = sqrt(4one(T)*T(pi)/(2*l+1))*rlm
	return rlm
end

rlm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable; norm::YLMNorm=Laplace{Float64}()) = Rlm{Float64}(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable; norm=norm)
