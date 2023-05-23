function s_cons()
    T = 1:Int(24 / TIME_STEP)
    T_size = size(T)[1] # number of periods per day

    load_base = (XLSX.readxlsx(XLSX_SUMMER_LOAD_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
    load_summer = copy(load_base)
    load_EV = (XLSX.readxlsx(XLSX_EV_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
    ELECTRIC_VEHICLES && (load_summer .+= load_EV)
    load_size = size(load_base)[2]
    divide_by = load_size ÷ T_size # integer divide
    averaged_load = [mean(load_summer[j, i:i+divide_by-1]) for j in 1:size(load_base)[1], i in 1:divide_by:T_size*divide_by]
    time_steps_lost = load_size - T_size * divide_by
    (time_steps_lost > 0) && (println("$time_steps_lost \"5min steps\" have been lost due to TIME_STEP value"))

    if WINTER
        load_winter = (XLSX.readxlsx(XLSX_WINTER_LOAD_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
        load_HP = (XLSX.readxlsx(XLSX_HP_PATH)[1][:])'[1:Nu_size, :] .* df_bus.S_D_mva[Nu]
        load_base = hcat(load_base, load_winter)
        ELECTRIC_VEHICLES && (load_winter .+= load_EV)
        HEAT_PUMPS && (load_winter .+= load_HP)
        averaged_load_winter = [mean(load_winter[j, i:i+divide_by-1]) for j in 1:size(load_winter)[1], i in 1:divide_by:T_size*divide_by]
        averaged_load = hcat(averaged_load, averaged_load_winter)
        T_size *= 2
    end

    scaling_factor = maximum(sum(load_base, dims=1))
    averaged_load *= PEAK_POWER / scaling_factor
    S_CONSUMPTION = [zeros(Ns_size, T_size); averaged_load[1:Nu_size, :]]
    return S_CONSUMPTION
end

include("parameters_key_BFM_1P.jl")
LAUNCH_SENSITIVITY_ANALYSIS = true
TIME_STEP = 1
(WINTER, ELECTRIC_VEHICLES, HEAT_PUMPS) = (false, false, false)
include("parameters_BFM_1P.jl")
LAUNCH_SENSITIVITY_ANALYSIS = false
# S_CONSUMPTION = s_cons()
# plot(sum(S_CONSUMPTION, dims=1)', label="summer")

ELECTRIC_VEHICLES = false
WINTER = true
S_CONSUMPTION = s_cons()
println(sum(S_CONSUMPTION, dims=1))
plot(sum(S_CONSUMPTION, dims=1)', label="summer+winter")

# ELECTRIC_VEHICLES = true
# S_CONSUMPTION = s_cons()
# plot!(sum(S_CONSUMPTION, dims=1)', label="summer+EV")

ELECTRIC_VEHICLES = true
WINTER = true
S_CONSUMPTION = s_cons()
println(maximum(sum(S_CONSUMPTION, dims=1)))
plot!(sum(S_CONSUMPTION, dims=1)', label="summer+winter+EV")

ELECTRIC_VEHICLES = false
WINTER = true
HEAT_PUMPS = true
S_CONSUMPTION = s_cons()
println(maximum(sum(S_CONSUMPTION, dims=1)))
plot!(sum(S_CONSUMPTION, dims=1)', label="summer+winter+HP")

ELECTRIC_VEHICLES = true
WINTER = true
HEAT_PUMPS = true
S_CONSUMPTION = s_cons()
println(maximum(sum(S_CONSUMPTION, dims=1)))
plot!(sum(S_CONSUMPTION, dims=1)', label="summer+winter+EV+HP")

# # hline!([PEAK_POWER], label="$(Int(PEAK_POWER)) MVA")
# # vline!([288.5], label="summer/winter", color=:gray)
# plot!(size=(1200, 400), legend=:topleft)
# xlabel!("Time steps of 5 minutes")
# ylabel!("Apparent power [MVA]")
# # annotate!(245, 20, "Summer day", :gray)
# # annotate!(2*288-250, 20, "Winter day", :gray)
# # change font of labels, hline & vline en pointillés, titres d'axes lisibles !
