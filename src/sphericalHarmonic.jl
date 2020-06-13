#normalization factor
abstract type YLMNorm{T} end

struct Schmidt{T} <: YLMNorm{T}
end

struct Laplace{T} <: YLMNorm{T}
end

struct Nonorm{T} <: YLMNorm{T}
end

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
function ylmKCoefficient(N::Nonorm{T}, l::Int64, m::Int64) where T
 return one(T)
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
function ylm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable; norm::YLMNorm{T}=Laplace{Float64}()) where T

  if abs(m) > l
    throw(DomainError(m,"-l <= m <= l expected, but m = $m and l = $l."))
  end

  p = (z^2 - 1)^l

  for i = 1:l+abs(m)
    c = i <= l ? 1/(2one(T)*i) : one(T)
    p = c*differentiate(p, z)
  end

  if m > 0
	  out=ylmKCoefficient(norm, l, m)*ylmCosSinPolynomial(m,x,y)*p
	  if norm!=Nonorm{T}()
		  out*=sqrt(2one(T))
	  end
    return out
  elseif m < 0
	  out=ylmKCoefficient(norm, l, abs(m))*ylmSinSinPolynomial(abs(m),x,y)*p
	  if norm!=Nonorm{T}()
		  out*=sqrt(2one(T))
	  end
    return out
  else
    return ylmKCoefficient(norm, l, 0)*p
  end
end

# multiplying r^l*ylm(x,y,z)
function rlylm(l::Int64, m::Int, x::Variable, y::Variable, z::Variable; norm::YLMNorm{T}=Laplace{Float64}()) where T
	p = ylm(l,m,x,y,z;norm=norm)
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
function rlm(l::Int64, m::Int64, x::Variable, y::Variable, z::Variable; norm::YLMNorm{T}=Laplace{Float64}()) where T
	rlm = rlylm(l,m,x,y,z; norm=norm)
	if norm!=Nonorm{T}()
		rlm = sqrt(4one(T)*T(pi)/(2*l+1))*rlm
	end
	return rlm
end
