module Solver
using DifferentialEquations, StaticArrays, JLD2, OffsetArrays, LinearAlgebra, TensorOperations, Dates

const tspan=(0.0,100.0)
const t_name=tspan[2]

const ⊕ = Base.FastMath.add_fast # fast versions for fma
const ⊖ = Base.FastMath.sub_fast
const ⊗ = Base.FastMath.mul_fast
const ⨸ = Base.FastMath.div_fast
const psi = 1.0
const g = 1.0
const k = sqrt(g*psi)
const p = (psi,g/2.0,k)
const N = 32
const N_eq = 2N+1


#=-------------------Combo Linear and quadratic term---------------------=#

let
odummy = OffsetArray(zeros(ComplexF64,N_eq,N_eq),-N:N,-N:N)
id = ones(ComplexF64, N_eq)
out = similar(id)
global Q_comb
function Q_comb(u,p::NTuple{3,Float64},t::Float64)
u = OffsetArray(u,-N:N)
ran = eachindex(u)
  for q in ran
    @simd for n in ran
        if abs(n-q) <= N
        @fastmath @inbounds odummy[q,n] = -(2⊗p[2]⊗p[1]⊗u[n]⊗conj(u[n-q]) ⊕
        (p[2])⊗conj(p[1])⊗u[n]⊗u[q-n])
          if q == n
           @fastmath @inbounds odummy[q,q] = (((abs(q)^2)⊗((p[3])^2⨸2)⊖p[2]⊗abs2(p[1]))⊗u[q]⊖p[2]⊗(p[1])^2⊗conj(u[-q])) - 2⊗p[2]⊗p[1]⊗u[n]⊗conj(u[n-q]) - (p[2])⊗conj(p[1])⊗u[n]⊗u[q-n]
          end
          # if q == 0
          #   @inbounds odummy[q,n] = 0.0 + 0.0im::Complex{Bool}
          # end
        end
    end
end
odummy1 = parent(odummy)
mul!(out,odummy1, id)
return out
end
end

#=-------------------Cubic term---------------------=#


let
odummy = OffsetArray(zeros(ComplexF64, N_eq,N_eq,N_eq),-N:N,-N:N,-N:N)
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
        @fastmath @inbounds odummy[q,n,m] = -(p[2])⊗u[n]⊗conj(u[m])⊗u[q-n+m]::ComplexF64
       end
      #= if q == 0
        @inbounds odummy[q,n,m] = 0.0 + 0.0im::Complex{Bool}
       end=#
     end
   end
end
odummy1 = parent(odummy)
@tensor out[q] = odummy1[q,n,m] * id[n,m]
return out
end
end


let
global result
function result(u,p::NTuple{3,Float64},t::Float64)
out = rdiv!(axpy!(1,Q_comb(u,p,t),C_mat(u,p,t)),im)
return SVector{N_eq,ComplexF64}(out)
end
end





function Diff()
global p
global N
global N_eq


# u1=[0.0+0.0*im for i in 1:N_eq]  # simple initial conditions
# u1[Int64((N_eq-1)/2)]=0.1+0.0*im
# u1[Int64((N_eq-1)/2)+1]=0.0+0.0*im
# u1[Int64((N_eq-1)/2)+2]=0.1+0.0*im
# for 2nd and 3rd bifurcation only values from -3 to 3 contribute
u1 = [  -0.0020564791627454377 + 0.003477416061042765im                           
    0.00817693525428735 + 0.012240759816549157im                           
   0.012603567910559872 + 0.10283443629478997im                            
    -0.5984905107747291 + 0.9135986296573508im                             
   0.012603567910637421 + 0.10283443629477225im                            
    0.00817693525429959 + 0.012240759816549963im                           
 -0.0020564791627472865 + 0.003477416061041004im ]
   for i in 1:29
    push!(u1, 0.0+0.0im)
    pushfirst!(u1, 0.0+0.0im)
  end
u0=SVector{N_eq}(u1)

Analytic_N=N


prob = ODEProblem(result,u0,tspan,p)
sol = solve(prob,RK4(),dt = 1.0*10^(-6),dense = false,saveat = 0.001)

path="C:/Users/visenkon/Desktop"
path1="C:/Users/visenkon/Desktop/JLD2/3nd bifurcation of Hamiltonian case with nodes -3 to 3, no forcing, alpha = 0, k= $(p[3]), g=$(2*p[2]), psi=$(p[1]), N=$Analytic_N, t_final = $(tspan[2])"
rm(path1,force=true,recursive=true)
mkdir(path1)
cd(path1)
@save "sol_tf=$(tspan[2]).jld2" sol
cd(path)
println(Dates.now())
return "success"
end
end
