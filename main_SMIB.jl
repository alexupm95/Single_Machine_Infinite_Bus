cd(dirname(@__FILE__));

#=
CODE FOR SOLVING THE SWING EQUATION FOR A SINGLE MACHINE CONNECTED TO AN INFINITE BUS

Author:      Alex Junior da Cunha Coelho
Supervisors: Luis Badesa Bernardo and Araceli Hernandez Bayo
Affiliation: Technical University of Madrid
August 2025
=#

# Load the necessary packages
using LinearAlgebra
using DifferentialEquations
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


# Time span
t_ini = 0.0    # Initial time
t_end = 2.0    # Final time
t_step = 0.001 # Time step to integrate the differential equations
tspan      = (t_ini, t_end)      # Tuple with the start time and end time for solving the differential equations
save_times = t_ini:t_step:t_end  # Save results every 0.01 seconds

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
    E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr, ω_s = p   # Input parameters
    δ, Δω = u                                                # δ and Δω
    du[1] = ω_s * Δω                                         # dδ/dt
    du[2] = (1 / (2*H)) * (P_m - P_e(δ, t, E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr) - D*(du[1]/ω_s))  # dΔω/dt
end

# Define and solve System or Differential Equations
u0   = [δ0, Δω0]                                                             # Vector of initial condition used to solve the differential equations
p    = (E, V, Vth, Xth, Xpostf, P_init, t_fault, t_cr, ω_s)                       # Tuple of input parameter used in the P_e function
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
fig_delta = plot!(xlims=[0.0, 1.5], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_delta)
savefig(fig_delta, "Delta.png")

# Plot P_ele(t)
fig_P_e = plot(time, Pet, label="P_ele", xlabel="Time (s)", ylabel="Power (p.u.)", title="Electrical Power vs Time", legend=:topleft, lw=3, lc=:orange)
fig_P_e = plot!([t_fault], seriestype=:vline, line=:dash, color=:red, label="t_f")
fig_P_e = plot!([t_fault + t_cr], seriestype=:vline, line=:dash, color=:black, label="t_cl")
fig_P_e = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_P_e = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_P_e = plot!(xlims=[0.0, 1.5], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_P_e)
savefig(fig_P_e, "Electrical Power.png")

# Plot ω(t) 
fig_ω = plot(time, ωt, xlabel="Time (s)", ylabel="ω (p.u.)", label="ω", title="Rotor Speed vs Time", legend=:topleft, lw=3, lc=:snow4)
fig_ω = plot!([t_fault], seriestype=:vline, line=:dash, color=:red, label="t_f")
fig_ω = plot!([t_fault + t_cr], seriestype=:vline, line=:dash, color=:black, label="t_cl")
fig_ω = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_ω = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_ω = plot!(xlims=[0.0, 1.5], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_ω)
savefig(fig_ω, "Rotor Speed.png")


# Plot f(t)
fig_f = plot(time, ft, xlabel="Time (s)", ylabel="f (Hz)", label="f", title="Frequency vs Time", legend=:topleft, lw=3, lc=:deepskyblue)
fig_f = plot!([t_fault], seriestype=:vline, line=:dash, color=:red, label="t_f")
fig_f = plot!([t_fault + t_cr], seriestype=:vline, line=:dash, color=:black, label="t_cl")
fig_f = plot!(size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15))
fig_f = plot!(fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=7mm)
fig_f = plot!(xlims=[0.0, 1.5], gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
display(fig_f)
savefig(fig_f, "Frequency.png")

println("=========================================================================================")
println("Solution saved and stored at: ", pwd())
println("=========================================================================================")



