module Solver
  using DifferentialEquations, JLD2, Dates
  include("linear_and_quadratic_subroutines.jl")
  include("paths.jl")
  
  function HamiltonianSolver(tspan::Tuple{Float64,Float64}, save_path::String = default_save_path, ϵ::ComplexF64 = 0.0 + 0.0im, modes_to_perturb::Array{Any} = [])
    save_solution_to = "$ϵ, $save_path" == "" ? save_path : "$main_path/JLD2/$ϵ, $save_path" 
    rm(save_solution_to, force = true, recursive = true)
    mkdir(save_solution_to) 
    cd(save_solution_to) 
    start_time = Dates.now()
    
    # Initial conditions
    
    initial=[0.0+0.0*im for i in 1:N_eq]
    initial[Int64((N_eq-1)/2)]=0.1+0.0*im
    initial[Int64((N_eq-1)/2)+1]=0.0+0.0*im
    initial[Int64((N_eq-1)/2)+2]=0.1+0.0*im
    
    # Perturb specified nodes by ϵ if ϵ and perturbed mode array are 
    # passed.
    if !isempty(modes_to_perturb) && !ϵ == 0.0 +0.0im
      for mode in modes_to_perturb
        initial[mode + N + 1] += ϵ
      end
    end
    
    initial_conditions = SVector{N_eq}(initial)
    
    # Solving the system using RK4, with specified initial conditions 
    # and step-size
    prob = ODEProblem(result,initial_conditions,tspan,p)
    solution = solve(prob,RK4(), dt = 1.0*10^(-6), dense = false, saveat = 0.001)
    
    # Saving final state of the system
    final_state = last(solution)
    open("final_state.txt","a") do io
      for state in final_state
        println(io,state)
      end
    end
    
    # Saving solution data for post-processing and displaying run-time
    @save "sol_tf=$(tspan[2]).jld2" solution
    cd(main_path)
    println("Run time: $(round((Dates.now() - start_time), Minute))")
    return final_state
  end
end
