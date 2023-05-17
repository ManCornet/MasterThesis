using OrderedCollections, Dates, PrettyTables, Latexify
include("parameters_key_BFM_1P.jl")
LAUNCH_SENSITIVITY_ANALYSIS = true

# /!\ MAKE SURE TO HAVE THE SAME ORDER IN THE LOOP BELOW
# /!\ PUT DEFAULT VALUES FIRST FOR EACH PARAMETER
parameters_values = OrderedDict(
    "WEIGHT_OP_LIMITS" => [1e-2, 1e-3, 1e-1],

    # # CONSUMPTION
    "WINTER" => [true], # true only
    "ELECTRIC_VEHICLES" => [false, true], # F/T
    "HEAT_PUMPS" => [false, true], # F/T

    # # PV
    "MAX_PV_CAPACITY_PER_NODE" => [0.4, 0.0, 0.8, 1.6],
    "PV_SCALE_SUMMER_WINTER" => [0.1],

    # # User related costs
    "IMP_ELECTRICITY_ENRG_COST" => [0.3, 0.6, 0.9],
    "IMP_ELECTRICITY_DSO_COST" => [0.1, 0.2, 0.3], # EXP_ELECTRICITY_DSO_COST =
    "EXP_ELECTRICITY_ENRG_COST" => [0.1, 0.2, 0.3],
    "GRID_CONNECTION_COST" => [80, 120, 160], # changements significatifs ?
)
param_names = [i for i in keys(parameters_values)]
values = [parameters_values[i] for i in param_names]
default_values = [i[1] for i in values]
EV_idx = [i for i in 1:length(parameters_values) if param_names[i] == "ELECTRIC_VEHICLES"][1]
HP_idx = [i for i in 1:length(parameters_values) if param_names[i] == "HEAT_PUMPS"][1]
MAX_PV_idx = [i for i in 1:length(parameters_values) if param_names[i] == "MAX_PV_CAPACITY_PER_NODE"][1]

value_table = []
results = []
k=0
# Iterators.product: goes through all possible combinations of the parameters above
for current_values in Iterators.product(values...)

    # only one parameter changes at a time (except EV + HP)
    current_vector = [i for i in current_values]
    diff_from_default = (current_vector .!= default_values)
    !((sum(diff_from_default) <= 1) ||
      ((sum(diff_from_default) == 2) && diff_from_default[EV_idx] && diff_from_default[HP_idx]) ||
      ((sum(diff_from_default) == 2) && diff_from_default[EV_idx] && diff_from_default[MAX_PV_idx]) ||
      ((sum(diff_from_default) == 2) && diff_from_default[HP_idx] && diff_from_default[MAX_PV_idx]) ||
      ((sum(diff_from_default) == 3) && diff_from_default[EV_idx] && diff_from_default[HP_idx] && diff_from_default[MAX_PV_idx])) &&
        continue

    (
        WEIGHT_OP_LIMITS,

        # CONSUMPTION
        WINTER,
        ELECTRIC_VEHICLES,
        HEAT_PUMPS,

        # PV
        MAX_PV_CAPACITY_PER_NODE,
        PV_SCALE_SUMMER_WINTER,

        # User related costs
        IMP_ELECTRICITY_ENRG_COST,
        IMP_ELECTRICITY_DSO_COST,
        EXP_ELECTRICITY_ENRG_COST,
        GRID_CONNECTION_COST,
        
        ) = current_values
    k+=1
    println(current_vector, " $k")

    include("parameters_BFM_1P.jl")
    param_table = Vector()
    push!(param_table, collect(param_names))
    push!(param_table, [i for i in current_values])
    append!(value_table, [param_table])

    try
        include("bilevel_BFM_1P.jl")
    catch
        append!(results, [["failure", ""]])
        continue
    else
        append!(results, [printed_tables])
    end

    pretty_table(current_vector', header=(param_names))
end

# INCLUDE results_tables.jl !!!
XLSX.openxlsx("sensitivity_analysis_" * (Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")) * ".xlsx", mode="w") do xf
    sheet_n = 1
    for (current_param, current_result) in zip(value_table, results)
        sheet_n >= 2 && XLSX.addsheet!(xf)
        sheet = xf[sheet_n]
        XLSX.writetable!(sheet, current_param, ["", ""])
        df = DataFrame(current_result[2:end, :], current_result[1, :])
        df = permutedims(df, 1)
        (ndims(current_result) > 1) && (df = df[!, [1, 3, 2]])
        XLSX.writetable!(sheet, df, anchor_cell=XLSX.CellRef("D1"))
        sheet_n = sheet_n + 1
    end
end


LAUNCH_SENSITIVITY_ANALYSIS = false;