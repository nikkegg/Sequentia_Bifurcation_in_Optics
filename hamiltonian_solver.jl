module Solver
  using DifferentialEquations, JLD2, Dates
  include("linear_and_quadratic_subroutines.jl")
  include("paths.jl")
  
  function HamiltonianSolver(save_path::String=default_save_path)
    save_solution_to = save_path == "" ? save_path : "$main_path/JLD2/$save_path" 
    rm(save_solution_to, force=true, recursive=true)
    mkdir(save_solution_to) 
    cd(save_solution_to) 
    global p
    global N
    global N_eq
    start_time = Dates.now()
    
    # Initial conditions
    
    initial=[0.0+0.0*im for i in 1:N_eq]
    initial[Int64((N_eq-1)/2)]=0.1+0.0*im
    initial[Int64((N_eq-1)/2)+1]=0.0+0.0*im
    initial[Int64((N_eq-1)/2)+2]=0.1+0.0*im
    
    # for i in 1:29
    #   push!(initial, 0.0+0.0im)
    #   pushfirst!(initial, 0.0+0.0im)
    # end
    
    initial_conditions=SVector{N_eq}(initial)
    Analytic_N=N
    prob = ODEProblem(result,initial_conditions,tspan,p)
    solution = solve(prob,RK4(), dt = 1.0*10^(-6), dense = false, saveat = 0.001)
    
    final_state = last(solution)
    open("final_state.txt","a") do io
      for state in final_state
        println(io,state)
      end
    end
    
    @save "sol_tf=$(tspan[2]).jld2" solution
    cd(main_path)
    println("Run time: $((Dates.now() - start_time) / 1000) seconds")
    return final_state
  end
end
