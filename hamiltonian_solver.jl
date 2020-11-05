module Solver
using DifferentialEquations, JLD2, Dates
include("linear_and_quadratic_subroutines.jl")

# Starting and ending simulation times
const tspan=(0.0,100.0)
const t_name=tspan[2]

# fast versions for fma
const ⊕ = Base.FastMath.add_fast 
const ⊖ = Base.FastMath.sub_fast
const ⊗ = Base.FastMath.mul_fast
const ⨸ = Base.FastMath.div_fast

  
  function Diff()
    global p
    global N
    global N_eq
  
  
    # u1=[0.0+0.0*im for i in 1:N_eq]
    # u1[Int64((N_eq-1)/2)]=0.1+0.0*im
    # u1[Int64((N_eq-1)/2)+1]=0.0+0.0*im
    # u1[Int64((N_eq-1)/2)+2]=0.1+0.0*im
    
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
