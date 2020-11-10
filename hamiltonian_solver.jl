module Solver
  using DifferentialEquations, JLD2, Dates
  include("linear_and_quadratic_subroutines.jl")
  include("paths.jl")
  
  # tspan - start and end simulation times
  # save_path - location for file to be saved at
  # ϵ - complex perturbation paramater if perturbation is needed
  # modes_to_perturb - array of modes to be perturbed. Normally, these
  # are selected as modes with highest amplitudes.
  function HamiltonianSolver(tspan::Tuple{Float64,Float64}, save_path::String = default_save_path, ϵ::ComplexF64 = 0.0 + 0.0im, modes_to_perturb = [], starting_values = [])
    save_solution_to = "epsilon = $ϵ, $save_path" == "" ? save_path : "$main_path/JLD2/epsilon = $ϵ, $save_path" 
    rm(save_solution_to, force = true, recursive = true)
    mkdir(save_solution_to) 
    cd(save_solution_to) 
    start_time = Dates.now()
    
    # Initial conditions be specified by hand
    if isempty(starting_values)
      initial = [0.0+0.0*im for i in 1:N_eq]
      initial[Int64((N_eq-1)/2)] = 0.1+0.0*im
      initial[Int64((N_eq-1)/2)+1] = 0.0+0.0*im
      initial[Int64((N_eq-1)/2)+2] = 0.1+0.0*im
    else
      initial = starting_values
    end
    
    # Perturb specified nodes by ϵ if ϵ and perturbed mode array is 
    # passed as an argument.
    if !isempty(modes_to_perturb) && !(ϵ == 0.0 + 0.0im)
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
    final_states = last(solution)
    open("final_state.txt","a") do io
      for state in final_states
        println(io,state)
      end
    end
    
    # Saving absolute values of the final state
    open("final_state_mag.txt","a") do io
      for state in final_states
        println(io,abs2(state))
      end
    end
    
    # Saving solution data for post-processing and displaying run-time
    @save "sol_tf=$(tspan[2]).jld2" solution
    cd(main_path)
    println("Run time: $(round((Dates.now() - start_time), Minute))")
    return final_states
  end
  
  function NthBifurcation(tspan::Tuple{Float64,Float64}, bifurcation_number::Int64, ϵ::ComplexF64, modes_to_perturb::Array{Int64}, starting_values, save_path::String = default_save_path)
    if isempty(starting_values)
      starting_values = HamiltonianSolver(tspan, save_path)
    end
    
    for iteration in 1:bifurcation_number
      solution = HamiltonianSolver(tspan, save_path * " bifurcation $(iteration + 1)", ϵ, modes_to_perturb, starting_values)
        starting_values = solution
      end
    return
  end

  function bifurcate_for(epsilon, tspan, save_path, modes_to_perturb, starting_values)
    for ϵ in epsilon
      HamiltonianSolver(tspan, save_path, ϵ, modes_to_perturb, starting_values)
    end
  end
end
