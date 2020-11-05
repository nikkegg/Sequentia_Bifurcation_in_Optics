module Solver
  using DifferentialEquations, JLD2, Dates
  include("linear_and_quadratic_subroutines.jl")
  include("paths.jl")
  
  function HamiltonianSolver(nodes_to_perturb::Array{Int64}, tspan::Tuple{Float64,Float64}, save_path::String = default_save_path, ϵ::ComplexF64 = 0.0 + 0.0im)
    save_solution_to = "$ϵ, $save_path" == "" ? save_path : "$main_path/JLD2/$ϵ, $save_path" 
    rm(save_solution_to, force=true, recursive=true)
    mkdir(save_solution_to) 
    cd(save_solution_to) 
    start_time = Dates.now()
    
    # Initial conditions
    
    initial=[0.0+0.0*im for i in 1:N_eq]
    initial[Int64((N_eq-1)/2)]=0.1+0.0*im
    initial[Int64((N_eq-1)/2)+1]=0.0+0.0*im
    initial[Int64((N_eq-1)/2)+2]=0.1+0.0*im
    
    # Perturb specified nodes by ϵ if ϵ and vector array are passed.
    if !isempty(nodes_to_perturb) && !ϵ == 0.0 +0.0im
      for node in nodes_to_perturb
        initial[nodes + N + 1] += ϵ
      end
    end
    
    initial_conditions=SVector{N_eq}(initial)
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
    println("Run time: $(round((Dates.now() - start_time, Minute))) seconds")
    return final_state
  end
end
