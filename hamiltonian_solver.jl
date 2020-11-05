module Solver
  using DifferentialEquations, JLD2, Dates
  include("linear_and_quadratic_subroutines.jl")

  # Starting and ending simulation times
  const tspan=(0.0,100.0)
  const t_name=tspan[2]
  
  function Diff()
    global p
    global N
    global N_eq
    start_time = Dates.now()
    
    # Initial conditions
    
    u1=[0.0+0.0*im for i in 1:N_eq]
    u1[Int64((N_eq-1)/2)]=0.1+0.0*im
    u1[Int64((N_eq-1)/2)+1]=0.0+0.0*im
    u1[Int64((N_eq-1)/2)+2]=0.1+0.0*im
    
    # for i in 1:29
    #   push!(u1, 0.0+0.0im)
    #   pushfirst!(u1, 0.0+0.0im)
    # end
    
    u0=SVector{N_eq}(u1)
    Analytic_N=N
    prob = ODEProblem(result,u0,tspan,p)
    sol = solve(prob,RK4(),dt = 1.0*10^(-6),dense = false,saveat = 0.001)
    path="C:/Users/visenkon/Desktop"
    path1="C:/Users/visenkon/Desktop/JLD2/testing refactored code"
    rm(path1,force=true,recursive=true)
    mkdir(path1)
    cd(path1)
    
    final_state = last(sol)
    open("final_state.txt","a") do io
      for state in final_state
        println(io,state)
      end
    end
    
    @save "sol_tf=$(tspan[2]).jld2" sol
    cd(path)
    println("Run time: $(Dates.now() - start_time)")
    return final_state
  end
end
