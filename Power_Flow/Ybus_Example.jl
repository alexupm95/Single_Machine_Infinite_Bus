cd(dirname(@__FILE__));

# Load the necessary packages
using LinearAlgebra, SparseArrays, DataFrames

#=
CODE FOR CALCULATING THE ADMITTANCE MATRIX (Ybus)

Author:      Alex Junior da Cunha Coelho
Affiliation: Technical University of Madrid
August 2025
=#


# Function to calculate the Ybus (admittance) Matrix
function Ybus_function(nbus::Int64, branches_df::DataFrame; kwargs...)
    Ybus = zeros(ComplexF64, nbus, nbus) # Instantiate the admittance matrix

    nbranches = size(branches_data)[1] # Number of branches in the system

    # Calculate the admittance matrix including the line data
    for i = 1:nbranches
        k = branches_df.from_bus[i] # Index from bus
        m = branches_df.to_bus[i]   # Index to bus
        if branches_df.status[i] == true

            display("($k, $m)")

            ykm = 1 / (branches_df.r[i] + 1im*branches_df.x[i])         # Series admittance
            bkm_sh = branches_df.bsh[i] / 2                             # Shunt admittance
            akm = branches_df.tap[i] == 0.0 ? 1.0 : branches_df.tap[i]  # Tap
            Φkm = deg2rad(branches_df.shift[i])                         # Shift in radians

            Ybus[k,k] = Ybus[k,k] + ((akm^2) * ykm) + (1im * bkm_sh)
            Ybus[k,m] = Ybus[k,m] - (akm * ykm * exp(-1im * Φkm))
            Ybus[m,k] = Ybus[m,k] - (akm * ykm * exp(1im * Φkm))
            Ybus[m,m] = Ybus[m,m] + ykm + (1im * bkm_sh)
        end
    end

    # Update the admittance matrix including the shunt elements at the buses, if applicable
    if haskey(kwargs, :df_bus)
        buses_df = kwargs[:df_bus]  # Getting the extra arguments, if applicable
        nbus_sh = size(buses_df)[1] # Number of buses with shunt elements

        for (i, id) in enumerate(buses_df.bus) 
            Ybus[id,id] = Ybus[id,id] + (1im * buses_df.bsh[i]) # Include the shunt components of the nodes
        end
    end

    return SparseArrays.sparse(Ybus)
end

#                      Lines and Transformers Data
# -----------------------------------------------------------------------------------------------------------------
# from_bus   to_bus     r (p.u.)       x (p.u.)      bsh_total (p.u.)    tap (p.u.)      shift (deg)        status
# ------------------------------------------------------------------------------------------------------------------
branches_data = [
   1	       4	          0	            0.0576	       0                    1.02            15           true
   2	       7	          0	            0.0625	       0                    0                0           true
   3	       9	          0	            0.0586	       0                    0                0           true
   4	       5	          0.01	        0.085	       0.176                0                0           true
   4	       6	          0.017	        0.092	       0.158                0                0           true
   5	       7	          0.032	        0.161	       0.306                0                0           true
   6	       9	          0.039	        0.17	       0.358                0                0           true
   7	       8	          0.0085	    0.072	       0.149                0                0           true
   8	       9	          0.0119	    0.1008	       0.209                0                0           true
]

#   Bus Data
# ------------
#  bus    bsh
# ------------
bus_data = [
    8     0.5
]

if Sys.iswindows()         # If system == Windows
    Base.run(`cmd /c cls`) # Clean the terminal
else                       # If system is based on Unix
    Base.run(`clear`)      # Clean the terminal
end

nbus = 9 # Number of active buses in the system

# Convert the Matrix into a Dataframe
df_branches = DataFrame(
    from_bus = Int.(branches_data[:,1]),
    to_bus   = Int.(branches_data[:,2]),
    r        = Float64.(branches_data[:,3]),
    x        = Float64.(branches_data[:,4]),
    bsh      = Float64.(branches_data[:,5]),
    tap      = Float64.(branches_data[:,6]),
    shift    = Float64.(branches_data[:,7]),
    status   = Bool.(branches_data[:,8])
)

# Convert the Matrix into a Dataframe
df_bus = DataFrame(
    bus      = Int.(bus_data[:,1]),
    bsh      = Float64.(bus_data[:,2]),
)


#= Call the function to calculate the admittance matrix 
NOTE: The matrix bus_data is optional, as the system may not have shunt reactive elements connected to the buses.
=#
Ybus =  Ybus_function(nbus, df_branches; df_bus) 
