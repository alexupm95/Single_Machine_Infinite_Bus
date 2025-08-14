cd(dirname(@__FILE__));

#=
CODE FOR SOLVING THE EQUAL AREA CRITERION AND CRITICAL CLEARING ANGLE AND TIME
FOR A SINGLE MACHINE CONNECTED TO AN INFINITE BUS

Author:      Alex Junior da Cunha Coelho
Supervisors: Luis Badesa Bernardo and Araceli Hernandez Bayo
Affiliation: Technical University of Madrid
August 2025
=#

# Load the necessary packages
using LinearAlgebra
using DifferentialEquations
using Trapz
using Plots
using Measures

cd(pwd()) # Load the current folder to save the figures

if Sys.iswindows() # If system == Windows
    Base.run(`cmd /c cls`) # Clean the terminal
else # If system is based on Unix
    Base.run(`clear`)     # Clean the terminal
end

# ================
# Input Parameters
# ================
H       = 3.0            # Inertia constant (s)
f       = 60             # Synchronous frequency (Hz)
ω_s     = 2π * f         # Synchronous speed (rad/s)
V       = 1.0            # Infinite bus voltage (p.u.)
Vth     = 0.3333         # Voltage of the Thévenin equivalent (p.u.)
E       = 1.1            # Internal EMF of the generator (p.u.), calculated from power balance
Xpref   = 0.1            # Equivalent reactance in Pre-Fault
Xth     = 0.09333        # Thévenin equivalent reactance during the Fault
Xpostf  = 0.12           # Equivalent reactance in Post-Fault
P_init  = 5.0            # Initial electrical power output (p.u.)
P_m     = P_init         # Mechanical power input (p.u.)
t_fault = 0.5            # Time o the fault
t_cr    = 0.12           # Time to clear the fault (seconds)
D       = 0.0            # Damping torque (p.u.)
δ0      = deg2rad(27.04) # Initial angle (δ₀)
Δω0     = 0.0            # Initial speed deviation (Δω₀)
δu = deg2rad(180.0 - 27.04) # Limit for unstable angle

# Calculate the maximum power in the pre-fault, during fault and post-fault stages
angles = collect(0:π/1000:π)  # Vector of angles
Pmax_pref  = (E * V) / Xpref  # Maximum electrical power output in the pre-fault period
Pmax_fault = (E * Vth) / Xth  # Maximum electrical power output in the fault period
Pmax_postf = (E * V) / Xpostf # Maximum electrical power output in the post-fault period

# Critical clearing angle (in radians)
δcct = acos((-P_m * δ0 + P_m * δu - Pmax_fault * cos(δ0) + Pmax_postf * cos(δu)) / (Pmax_postf - Pmax_fault))

println("==============================")
println("Critical Angle: "*string(round(rad2deg(δcct), digits=3))*" degrees")
println("==============================\n")

# =======================================================
# Calculate the Equal Area Criterion to confirm stability
# =======================================================
indices_fault   = findall(x -> δ0 ≤ x ≤ δcct, angles)       # Indices of angles in the fault interval
δ_fault         = angles[indices_fault]                     # Angles in the fault interval
Pacc_A1         = P_m .* ones(Float64, length(indices_fault)) .- Pmax_fault .* sin.(δ_fault) # Accelerating power Area 1
# integral_Area_1 = Trapz.trapz(δ_fault, Pacc_A1)                       # Integral of Area 1 during the fault (using Trapezoidal integration)
integral_Area_1 = P_m * (δcct - δ0) - Pmax_fault*(-cos(δcct) + cos(δ0)) # Integral of Area 1 during the fault (using direct formula)


indices_postf = findall(x -> δcct ≤ x ≤ δu, angles)       # Indices of angles in the post-fault interval
δ_postf       = angles[indices_postf]                     # Angles in the post-fault interval
Pacc_A2       = -P_m .* ones(Float64, length(indices_postf)) .+ Pmax_postf .* sin.(δ_postf)  # Accelerating power Area 1
# integral_Area_2 = Trapz.trapz(δ_postf, Pacc_A2)                         # Integral of Area 2 during the fault (using Trapezoidal integration)
integral_Area_2 = - P_m * (δu - δcct) + Pmax_postf*(-cos(δu) + cos(δcct)) # Integral of Area 2 during the fault (using direct formula)

# Print the values of the integrals
println("==========================")
println("Integral in Area 1: "*string(round(integral_Area_1, digits=4))*"")
println("Integral in Area 2: "*string(round(integral_Area_2, digits=4))*"")
println("==========================\n")


# Plot the curves of power output
fig_power = plot(rad2deg.(angles), Pmax_pref .* sin.(angles), label="P_pref", legend=:topright, lw=2, lc=:limegreen)
fig_power = plot!(rad2deg.(angles), Pmax_fault .* sin.(angles), label="P_fault", lw=2, lc=:orange)
fig_power = plot!(rad2deg.(angles), Pmax_postf .* sin.(angles), label="P_postf", lw=2, lc=:deepskyblue)
fig_power = plot!([P_m], seriestype=:hline, line=:dash, color=:red, label="P_mec")
fig_P_ele = plot!(rad2deg.(δ_fault), Pmax_fault .* sin.(δ_fault), fillrange=P_m, fillalpha=0.3, color=:lime, label="", legend=:topright)
fig_P_ele = plot!(rad2deg.(δ_postf), Pmax_postf .* sin.(δ_postf), fillrange=P_m, fillalpha=0.3, color=:magenta, label="", legend=:topright)
fig_power = plot!([rad2deg(δ0), rad2deg(δ0)], [0, P_m],label="δ₀", lw=2, ls=:dash, lc=:black)
fig_power = plot!([rad2deg(δcct), rad2deg(δcct)], [Pmax_fault*sin(δcct), Pmax_postf*sin(δcct)],label="δcct", lw=2, ls=:dash, lc=:snow4)
fig_power = plot!([rad2deg(δu), rad2deg(δu)], [0, P_m],label="δᵤ", lw=2, ls=:dash, lc=:red)
fig_power = annotate!([(55, 4.1, ("Area 1", 22, :black, :center, "courier")), (115, 6.3, ("Area 2", 22, :black, :center,"courier"))])
fig_power = plot!(xlabel="Angle (deg)", ylabel="Power (p.u.)", title="Power-Angle Curves")
fig_power = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_power = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_power = plot!(xlims=[0.00, 190], xticks = [30, 60, 90, 120, 150, 180], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_power)
savefig(fig_power, "EAC.png")

# =============================================================
# Solve the Swing Equation to obtain the Critical Clearing Time
# =============================================================

# Time span
t_ini      = 0.0                 # Initial time
t_end      = 2.0                 # Final time
t_step     = 0.001               # Time step to integrate the differential equations
tspan      = (t_ini, t_end)      # Tuple with the start time and end time for solving the differential equations
save_times = t_ini:t_step:t_end  # Save results every 0.01 seconds

# Callback to monitor the Critical Clearing Time
function Callback_Monitor_CCT(u, t, integrator)
    E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr, ω_s, δcct = integrator.p
    δ = u[1]

    if (t ≥ t_fault) && (δcct - δ) ≤ 0.001 # If δcct - δ is lower than the tolerance, the program stops to solve
        return true
    
    else
        return false
    end
end

# Function to calculate the electrical power depending on the values of δ and t
function P_e(δ, t, E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr)
    if 0.0 ≤ t < t_fault # Pre-fault
        return P_init

    elseif t_fault ≤ t < (t_fault + t_cr)
        return (Vth * E * sin(δ)) / Xth  # During fault - bus 3 is short-circuited

    else
        return (V * E * sin(δ)) / Xpostf # Post-fault
    end
end

# Swing equation
function swing_eq!(du, u, p, t)
    E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr, ω_s, δcct = p  # Input parameters
    δ, Δω = u                                                     # δ and Δω
    du[1] = ω_s * Δω                                              # dδ/dt
    du[2] = (1 / (2*H)) * (P_m - P_e(δ, t, E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr) - D*(du[1]/ω_s))  # dΔω/dt
end

# Define and solve System or Differential Equations
u0   = [δ0, Δω0]                                                             # Vector of initial condition used to solve the differential equations
p    = (E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr, ω_s, δcct)                            # Tuple of input parameter used in the P_e function

affect!(integrator) = terminate!(integrator)                          # If condition is met, terminate the integrator
cb_cct              = DiscreteCallback(Callback_Monitor_CCT, affect!) # Callback -> Terminate the ODE solver when the condition is met         
prob                = ODEProblem(swing_eq!, u0, tspan, p)                                                    # Define the problem
sol_cct             = solve(prob, Trapezoid(), callback=cb_cct, reltol=1e-6, abstol=1e-8, saveat=save_times) # Solve the system

cct = sol_cct.t[end] # Get the Critical Clearing Time (last position of the time solution of the Swing Equation)

# Print the values of the integrals
println("=================================")
println("Critical Clearing Time: "*string(round(cct, digits=3))*" sec")
println("=================================\n")

# ==========================================================
# Solve the Swing Equation to obtain the curves for plotting
# ==========================================================
t_cr =  round(cct; digits=3) - t_fault # Clearing time = Critical Clearing Time

# Define and solve System or Differential Equations
u0   = [δ0, Δω0]                                                             # Vector of initial condition used to solve the differential equations
p    = (E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr, ω_s, δcct)            # Tuple of input parameter used in the P_e function
prob = ODEProblem(swing_eq!, u0, tspan, p)                                   # Define the problem
sol  = solve(prob, Trapezoid(), reltol=1e-6, abstol=1e-8, saveat=save_times) # Solve the system

# Save the solution in vectors
time = sol.t                    # Vector of time (sec)
δt   = rad2deg.(sol[1, :])      # Vector of delta over time (deg)
ωt  = 1.0 .+ sol[2, :]          # Vector of rotor speed over time (p.u.)
ft  = f .* ωt                   # Electrical frequency (Hz)
Pet  = @. P_e(sol[1, :], sol.t, E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr) # Vector of electrical power over time (p.u.)

# =========================
# Plot the output variables
# =========================

# Plot δ(t)
fig_delta = plot(time, δt, label="δ", xlabel="Time (s)", ylabel="δ (deg)", title="δ vs Time", legend=:topleft, lw=3, lc=:limegreen)
fig_delta = plot!([t_fault], seriestype=:vline, line=:dash, color=:red, label="t_f", lw=2.0)
fig_delta = plot!([t_fault + t_cr], seriestype=:vline, line=:dash, color=:black, label="t_cl", lw=2.0)
fig_delta = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_delta = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_delta = plot!(xlims=[0.0, 2.0], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_delta)
savefig(fig_delta, "Delta.png")

# Plot P_ele(t)
fig_P_e = plot(time, Pet, label="P_ele", xlabel="Time (s)", ylabel="Power (p.u.)", title="Electrical Power vs Time", legend=:topleft, lw=3, lc=:orange)
fig_P_e = plot!([t_fault], seriestype=:vline, line=:dash, color=:red, label="t_f")
fig_P_e = plot!([t_fault + t_cr], seriestype=:vline, line=:dash, color=:black, label="t_cl")
fig_P_e = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_P_e = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_P_e = plot!(xlims=[0.0, 2.0], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_P_e)
savefig(fig_P_e, "Electrical Power.png")

# Plot ω(t) 
fig_ω = plot(time, ωt, xlabel="Time (s)", ylabel="ω (p.u.)", label="ω", title="Rotor Speed vs Time", legend=:topleft, lw=3, lc=:snow4)
fig_ω = plot!([t_fault], seriestype=:vline, line=:dash, color=:red, label="t_f")
fig_ω = plot!([t_fault + t_cr], seriestype=:vline, line=:dash, color=:black, label="t_cl")
fig_ω = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_ω = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_ω = plot!(xlims=[0.0, 2.0], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_ω)
savefig(fig_ω, "Rotor Speed.png")


# Plot f(t)
fig_f = plot(time, ft, xlabel="Time (s)", ylabel="f (Hz)", label="f", title="Frequency vs Time", legend=:topleft, lw=3, lc=:deepskyblue)
fig_f = plot!([t_fault], seriestype=:vline, line=:dash, color=:red, label="t_f")
fig_f = plot!([t_fault + t_cr], seriestype=:vline, line=:dash, color=:black, label="t_cl")
fig_f = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_f = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_f = plot!(xlims=[0.0, 2.0], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_f)
savefig(fig_f, "Frequency.png")

println("=========================================================================================")
println(" Figures stored at: ", pwd())
println("=========================================================================================")
