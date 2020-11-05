using OffsetArrays, TensorOperations, StaticArrays, LinearAlgebra


# psi - strength of initial Einstein-Bose condensate being perturbed
# g non-linearity (Kerr) parameter
# k wavenumber of the most unstable mode
# N - number of modes in Fourier expansion
# N_eq = total number of equation, as Fourier series goes from -N to N

const psi = 1.0
const g = 1.0
const k = sqrt(g * psi)
const p = (psi, g/2.0, k)
const N = 32
const N_eq = 2N + 1

let
  coefficient_array = OffsetArray(zeros(ComplexF64, N_eq, N_eq), -N:N,   -N:N)
  id = ones(ComplexF64, N_eq)
  out = similar(id)
  global Q_comb
  function Q_comb(u, p::NTuple{3,Float64}, t::Float64)
    u = OffsetArray(u,-N:N)
    ran = eachindex(u)
    for q in ran
      @simd for n in ran
        if abs(n-q) <= N
          @fastmath @inbounds coefficient_array[q,n] = -(2⊗p[2]⊗p[1]⊗u[n]⊗conj(u[n-q]) ⊕ (p[2])⊗conj(p[1])⊗u[n]⊗u[q-n])
          if q == n
            @fastmath @inbounds coefficient_array[q,q] = (((abs(q)^2)⊗(  (p[3])^2⨸2)⊖p[2]⊗abs2(p[1]))⊗u[q]⊖p[2]⊗(p[1])^2⊗conj(u[-q])) - 2⊗p[2]⊗p[1]⊗u[n]⊗conj(u[n-q]) - (p[2])⊗conj(p[1])⊗u[n]⊗u[q-n]
          end
        # if q == 0
          # @inbounds coefficient_array[q,n] = 0.0 +     0.0im::Complex{Bool}
        # end
        end
      end
    end
    coefficients = parent(coefficient_array)
    mul!(out,coefficients, id)
    return out
  end
end
  
  #=-------------------Cubic term---------------------=#
  
  
let
  coefficient_array = OffsetArray(zeros(ComplexF64, N_eq,N_eq,N_eq),-N:N,-N:N,-N:N)
  id = ones(ComplexF64, N_eq,N_eq)
  out = zeros(ComplexF64,N_eq)
  global C_mat
  function C_mat(u,p::NTuple{3,Float64},t::Float64)
    u = OffsetArray(u,-N:N)
    ran = eachindex(u)
    for q in ran
      for n in ran
        @simd for m in ran
          if abs(q-n+m) <= abs(N)
            @fastmath @inbounds coefficient_array[q,n,m] = -(p[2])⊗u[n]⊗conj(u[m])⊗u[q-n+m]::ComplexF64
          end
       #= if q == 0
            @inbounds coefficient_array[q,n,m] = 0.0 + 0.0im::Complex{Bool}
          end=#
        end
      end
    end
    coefficients = parent(coefficient_array)
    @tensor out[q] = coefficients[q,n,m] * id[n,m]
    return out
  end
end
  
  
let
  global result
  function result(u, p::NTuple{3,Float64}, t::Float64)
    out = rdiv!(axpy!(1, Q_comb(u,p,t), C_mat(u,p,t)), im)
    return SVector{N_eq, ComplexF64}(out)
  end
end
