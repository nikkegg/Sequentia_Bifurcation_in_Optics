include("C:/Users/visenkon/Desktop/hamiltonian_solver.jl")          
module Plt
using Plots
pyplot()
import PyPlot
function plots(N::Int64,psi::Float64,g::Float64,spec_plot::Array{T},α::Float64,ν::Float64, to_plot::Bool) where T<:Int64 
path1 = "C:/Users/visenkon/Desktop/Corrected Phd write-up figure 1, alph = 0, k= 1.0, g=1.0, psi=1.0, N=32, t_final = 100.0"
cd(path1)
tspan=(0.0, 100.0)
JLD2.@load "sol_tf=$(tspan[2]).jld2" sol
if to_plot
t_name=tspan[2]
Analytic_N=Int64(floor((N-1)/2))
path="C:/Users/visenkon/"
# path1="C:/Users/visenkon/Desktop/JLD2/Dissipative solution for alpha = , g=$g, psi=$psi, N=$Analytic_N, t_final = $(tspan[2])"#"C:/Users/visenkon/Desktop/JLD2 tfinal =400 t0 = 0.1,for g=$g, psi=$psi, N=$Analytic_N,t=$t_name"
realp=[real.(sol[i,:]) for i in 1:N]
imp=[imag.(sol[i,:]) for i in 1:N]
out_mag=[abs2.(sol[i,:]) for i in 1:N]
out_t=sol.t[:]
realp=[real.(sol[i,:]) for i in 1:N]
imp=[imag.(sol[i,:]) for i in 1:N]
sol_cons = [setindex(sol[i],sol[i][Analytic_N + 1] + psi, Analytic_N + 1) for i in 1:length(out_t)]
sol_mag = [abs2.(sol_cons[i]) for i in 1:length(out_t)]   
for_phase = realp
for_phase[33][1:end] .+= 1.0
path2=string("C:/Users/visenkon/Desktop/Plotsss/Hyper, 3rd, nu = $ν ")
mkdir(path2)
cd(path2)
#= 
magnitutd eplots on different x;lims **************************************************************************
=#
kmin=min(spec_plot...)
magnitude1=plot(out_t,out_mag[spec_plot[1]], xlims = (0, 50), title = L"Amplitudes of $A_k$", #=label = "k=$(spec_plot[1] - Analytic_N -1)"*L"k^*",=# 
legend = false,#=:outertopright=# xaxis="time, t",#=legendfontsize=15, legendfont = "serif", foreground_color_legend = nothing, background_color_legend = nothing,=# titlefont = "serif", 
titlefontsize = 20, yaxis=L"|A_k|^2",lw=1.3, grid = false,  xguidefontsize=16, yguidefontsize=16, xguidefontfamily = "serif", yguidefontfamily = "serif")
for k in spec_plot[2:end]
plot!(out_t,out_mag[k],label="k=$(k - Analytic_N - 1)"*L"k^*",lw=1.3)
end


magnitude2=plot(out_t,out_mag[spec_plot[1]], xlims = (0, 10), title = L"Amplitudes of $A_k$", label = "k=$(spec_plot[1] - Analytic_N -1)"*L"k^*", 
legend = :outertopright, xaxis="time, t",legendfontsize=15, legendfont = "serif", titlefont = "serif", 
titlefontsize = 20, yaxis=L"|A_k|^2",lw=1.3, grid = false, background_color_legend = nothing, xguidefont=font(16), yguidefont=font(16))
for k in spec_plot[2:end]
plot!(out_t,out_mag[k],label="k=$(k - Analytic_N - 1)"*L"k^*",lw=1.3)
end

magnitude3=plot(out_t,out_mag[spec_plot[1]], title = L"Amplitudes of $A_k$", label = "k=$(spec_plot[1] - Analytic_N -1)"*L"k^*", 
legend = false, xaxis="time, t",legendfontsize=15, legendfont = "serif",titlefont = "serif", 
foreground_color_legend = nothing, background_color_legend = nothing,
titlefontsize = 20, yaxis=L"|A_k|^2",lw=1.3, grid = false,  xguidefont=font(16), yguidefont=font(16))
for k in spec_plot[2:end]
plot!(out_t,out_mag[k],label="k=$(k - Analytic_N - 1)"*L"k^*",lw=1.3)
end

#=
plots of imaginary parts *********************************************************************
=#

for k in spec_plot
plot(out_t,imp[k][1:end], title="k=$(k - Analytic_N - 1)"*L"$k^*$", 
	xlabel = "time, t",ylabel= L"Im($A_k$)", lw=1.0,
	 legend = false, xaxis="time, t", titlefont = "serif",
    titlefontsize = 20, grid = false,  xguidefont=font(16), yguidefont=font(16), linecolor = :orange)
savefig("imaginary_part_node$(k - Analytic_N - 1)k.eps")
end

# #=
# Phase plots  *****************************************************************************8
# =#

for k in spec_plot
plot(for_phase[k][1:8:end],imp[k][1:8:end], title="k=$(k - Analytic_N - 1)"*L"$k^*$", xlabel = L"Re($A_k$)",
	ylabel= L"Im($A_k$)", lw=1.0, grid=false, leg=false, linecolor = :orange, titlefont = "serif",
    titlefontsize = 20, xguidefont=font(16), yguidefont=font(16))
savefig("phase_plot_node_$(k - Analytic_N - 1)k.eps")
end

# #=
# Wave action plots *********************************************************************************
# =#

c_plot = sum.(sol_mag[:])
pcons=plot(out_t, c_plot,lw=2.0,grid=false,title="Wave action",
	xlabel="time, t",ylabel=L"$\sum_{i=-N}^{N} |A_i|^2$",
	leg=false,titlefont = "serif",titlefontsize = 20, xguidefont=font(20), 
	yguidefont=font(20), linecolor = :green)

#=
Spectrum plots ***************************************************************************************
=#

# /

t20=[out_mag[:,1][:][i][1] for i in 1:N]
t40=[out_mag[:,1][:][i][200000] for i in 1:N]
t60=[out_mag[:,1][:][i][300000] for i in 1:N]
t80=[out_mag[:,1][:][i][400000] for i in 1:N]
t100=[out_mag[:,1][:][i][500000] for i in 1:N]
t200=[out_mag[:,1][:][i][600000] for i in 1:N]
t300=[out_mag[:,1][:][i][700000] for i in 1:N]
t400=[out_mag[:,1][:][i][800000] for i in 1:N]
# t500=[out_mag[:,1][:][i][900000] for i in 1:N]
# t600=[out_mag[:,1][:][i][1000000] for i in 1:N]
# t700=[out_mag[:,1][:][i][1100000] for i in 1:N]
# t800=[out_mag[:,1][:][i][1200000] for i in 1:N]
# t900=[out_mag[:,1][:][i][1300000] for i in 1:N]
# t1000=[out_mag[:,1][:][i][1400000] for i in 1:N]
# t1100=[out_mag[:,1][:][i][1500000] for i in 1:N]
# t1200=[out_mag[:,1][:][i][1600000] for i in 1:N]
offset_t20=OffsetArray(t20,-Analytic_N:Analytic_N)
offset_t40=OffsetArray(t40,-Analytic_N:Analytic_N)
offset_t60=OffsetArray(t60,-Analytic_N:Analytic_N)
offset_t80=OffsetArray(t80,-Analytic_N:Analytic_N)
offset_t100=OffsetArray(t100,-Analytic_N:Analytic_N)
offset_t200=OffsetArray(t200,-Analytic_N:Analytic_N)
offset_t300=OffsetArray(t300,-Analytic_N:Analytic_N)
offset_t400=OffsetArray(t400,-Analytic_N:Analytic_N)
# offset_t500=OffsetArray(t40,-Analytic_N:Analytic_N)
# offset_t600=OffsetArray(t60,-Analytic_N:Analytic_N)
# offset_t700=OffsetArray(t80,-Analytic_N:Analytic_N)
# offset_t800=OffsetArray(t100,-Analytic_N:Analytic_N)
# offset_t900=OffsetArray(t200,-Analytic_N:Analytic_N)
# offset_t1000=OffsetArray(t300,-Analytic_N:Analytic_N)
# offset_t1100=OffsetArray(t400,-Analytic_N:Analytic_N)
# offset_t1200=OffsetArray(t400,-Analytic_N:Analytic_N)
offset_spec=[offset_t20,offset_t40,offset_t60,offset_t80,offset_t100,
offset_t200,offset_t300,offset_t400]
#,offset_t500,offset_t600,offset_t700,offset_t800,offset_t900,offset_t1000,
# offset_t1100,offset_t1200]
# offset_spec=[offset_t20,offset_t40,offset_t80,offset_t200,offset_t400]
l_offset_spec=length(offset_spec)

spec20=Vector{Float64}(undef,Analytic_N)
spec40=Vector{Float64}(undef,Analytic_N)
spec60=Vector{Float64}(undef,Analytic_N)
spec80=Vector{Float64}(undef,Analytic_N)
spec100=Vector{Float64}(undef,Analytic_N)
spec200=Vector{Float64}(undef,Analytic_N)
spec300=Vector{Float64}(undef,Analytic_N)
spec400=Vector{Float64}(undef,Analytic_N)
# spec500=Vector{Float64}(undef,Analytic_N)
# spec600=Vector{Float64}(undef,Analytic_N)
# spec700=Vector{Float64}(undef,Analytic_N)
# spec800=Vector{Float64}(undef,Analytic_N)
# spec900=Vector{Float64}(undef,Analytic_N)
# spec1000=Vector{Float64}(undef,Analytic_N)
# spec1100=Vector{Float64}(undef,Analytic_N)
# spec1200=Vector{Float64}(undef,Analytic_N)

spec=[spec20,spec40,spec60,spec80,spec100,spec200,spec300,spec400]#,spec500,spec600,spec700,
#spec800,spec900,spec1000,spec1100,spec1200]
# # spec=[spec20,spec40,spec80,spec200,spec400]
spec_label=["t=0","t=200","t=300","t=400","t=500","t=600","t=700","t=800"]#"t=900","t=1000","t=1100",
#"t=1200","t=1300","t=1400","t=1500","t=1600"]
# spec_label=["t=0","t=20","t=40","t=60","t=80"]
for k=1:length(spec)
for i=1:Analytic_N
spec[k][i]=offset_spec[k][-i]+offset_spec[k][i]
end
end
lc = [:red, :green, :blue, :orange, :purple, :black, :pink, :cyan] 
mod_n=[i for i in 1:Analytic_N]
spectrum=plot(mod_n,spec[1], title = "Spectrum", legend = :outertopright, legendfontsize=15, 
	legendfont = "serif",foreground_color_legend = nothing, 
	background_color_legend = nothing, titlefont = "serif", titlefontsize = 20, 
	grid=false,xlabel="wavenumber, k", ylabel=L"n_k =|A_k|^2 + |A_{-k}|^2",
	xaxis=:log10, yaxis=:log10, label="t=0", lw=2,xguidefont=font(20), yguidefont=font(20), 
	linecolor = lc[1], ylims = (10^-12,10^4), xlims = (10^0.0, 10^1.4))
for i=2:length(spec)
plot!(mod_n,spec[i],lw=1.5,xaxis=:log10, yaxis=:log10, label=spec_label[i],grid=false,
	xguidefont=font(16), yguidefont=font(16), linecolor = lc[i])
end



savefig(pcons, "conservation.eps")
savefig(magnitude1,"mag1.eps")
# savefig(magnitude2,"mag2.eps")
savefig(magnitude3,"mag3.eps")
savefig(spectrum,"spectrum.eps")
cd(path)
else
return sol
end
end

nu = [0.0027,0.0028,0.0029,0.003, 0.0031]
function multplot()
for ν in nu
	plots(41, 1.0, 1.0, collect(1:41), 0.1, ν, true)
end
end
end
