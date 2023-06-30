using XLSX, DataFrames, Plots, Statistics
# include("parameters.jl")

file_name = "output_2023-04-27_17-31-35"
file_path = splitdir(pwd())[1] * "\\bilevel\\" * file_name * ".xlsx"
file = XLSX.readxlsx(file_path)


# variables import
S_CONSUMPTION = Matrix{Float64}(XLSX.getdata(file["S_CONSUMPTION"])[2:end, 2:end])
P_CONSUMPTION = Matrix{Float64}(XLSX.getdata(file["P_CONSUMPTION"])[2:end, 2:end])
Q_CONSUMPTION = Matrix{Float64}(XLSX.getdata(file["Q_CONSUMPTION"])[2:end, 2:end])
S_MAX_LINE = Matrix{Float64}(XLSX.getdata(file["S_MAX_LINE"])[2:end, 2:end])
Alpha = Matrix{Int}(round.(XLSX.getdata(file["Alpha"])[2:end, 2:end]))
P_line = Matrix{Float64}(XLSX.getdata(file["P_line"])[2:end, 2:end])
Q_line = Matrix{Float64}(XLSX.getdata(file["Q_line"])[2:end, 2:end])
p_pv = Matrix{Float64}(XLSX.getdata(file["p_pv"])[2:end, 2:end])
q_pv = Matrix{Float64}(XLSX.getdata(file["q_pv"])[2:end, 2:end])
S_sub = round.(Vector{Float64}(XLSX.getdata(file["S_sub"])[2:end, 2]), digits=3)
user_costs = Matrix{Float64}(XLSX.getdata(file["user_costs"])[2:end, 2:end])
s_grid_max = Matrix{Float64}(XLSX.getdata(file["s_grid_max"])[2:end, 2:end])
DSO_costs = Int(round(XLSX.getdata(file["DSO_costs"])[2]))
P_line = Matrix{Float64}(XLSX.getdata(file["P_line"])[2:end, 2:end])
Q_line = Matrix{Float64}(XLSX.getdata(file["Q_line"])[2:end, 2:end])
S_MAX_build_line = Vector{Float64}(XLSX.getdata(file["S_MAX_build_line"])[2:end, 2])
p_imp = Matrix{Float64}(XLSX.getdata(file["p_imp"])[2:end, 2:end])



# sets reconstruction
L_size = size(P_line)[1]
L = 1:L_size
T_size = size(P_line)[2]
T = 1:T_size
Ns_size = size(S_sub)[1]
Ns = 1:Ns_size
Nu_size = size(p_pv)[1]
Nu = (1:Nu_size) .+ Ns_size


# KEY INDICATORS

# self sufficiency
p_self_sufficiency = min.(p_pv, P_CONSUMPTION[Nu, :]) ./ P_CONSUMPTION[Nu, :]
q_self_sufficiency = ifelse.(q_pv .< Q_CONSUMPTION[Nu, :], q_pv, Q_CONSUMPTION[Nu, :]) ./ Q_CONSUMPTION[Nu, :]
q_self_sufficiency = ifelse.(q_self_sufficiency .< 0, 0, q_self_sufficiency)

b_range = [-0.01, 0.01, 0.25, 0.5, 0.75, 0.99, 1.01]

hist1 = histogram(
    reshape(p_self_sufficiency, :, 1),
    bins=b_range,
    normalize=:probability,
    legend=false)
xlabel!("p_self_sufficiency /1")
ylabel!("Frequency /1")
# ylims!(0, 0.8)
vline!([mean(p_self_sufficiency)])
display(hist1)

hist2 = histogram(
    reshape(q_self_sufficiency, :, 1),
    bins=b_range,
    normalize=:probability,
    legend=false)
xlabel!("q_self_sufficiency /1")
ylabel!("Frequency /1")
# ylims!(0, 0.8)
vline!([mean(q_self_sufficiency)])
display(hist2)

# self consumption
p_self_consumption = min.(p_pv, P_CONSUMPTION[Nu, :]) ./ p_pv
q_self_consumption = ifelse.(q_pv .< Q_CONSUMPTION[Nu, :], q_pv, Q_CONSUMPTION[Nu, :]) ./ q_pv
q_self_consumption = ifelse.(q_self_consumption .< 0, 0, q_self_consumption)

hist3 = histogram(
    reshape(p_self_consumption, :, 1),
    bins=b_range,
    normalize=:probability,
    legend=false)
xlabel!("p_self_consumption /1")
ylabel!("Frequency /1")
ylims!(0, 1)
vline!([mean(p_self_consumption)])
display(hist3)

hist4 = histogram(
    reshape(q_self_consumption, :, 1),
    bins=b_range,
    normalize=:probability,
    legend=false)
xlabel!("q_self_consumption /1")
ylabel!("Frequency /1")
ylims!(0, 1)
vline!([mean(q_self_consumption)])
display(hist4)


# user costs
hist5 = histogram(
    reshape(user_costs, :, 1),
    bins=5,
    normalize=:probability,
    legend=false)
xlabel!("user costs [EUR]")
ylabel!("Frequency /1")
# ylims!(0,0.8)
vline!([mean(user_costs)])
display(hist5)


# s grid max
hist6 = histogram(
    reshape(s_grid_max, :, 1),
    bins=5,
    normalize=:probability,
    legend=false)
xlabel!("user grid capacity [kVA]")
ylabel!("Frequency /1")
# ylims!(0,0.8)
vline!([mean(s_grid_max)])
display(hist6)


# p imp
p_imp_user = mean(p_imp, dims=2)
hist7 = histogram(
    reshape(p_imp_user, :, 1),
    # bins=5,
    normalize=:probability,
    legend=false)
xlabel!("p imported, mean over a day [kW]")
ylabel!("Frequency /1")
# ylims!(0,0.8)
vline!([mean(p_imp_user)])
display(hist7)


# line load
S_line = sqrt.(P_line .^ 2 + Q_line .^ 2)
line_load = S_line ./ S_MAX_LINE
line_load = line_load[vec(mapslices(col -> any(col .!= 0), line_load, dims=2)), :] # removes lines full of zeros
hist8 = histogram(
    reshape(line_load, :, 1),
    bins=6,
    normalize=:probability,
    legend=false)
xlabel!("line load [kVA]")
ylabel!("Frequency /1")
# ylims!(0,0.8)
vline!([mean(line_load)])
display(hist8)


println(file_name)
println("DSO_costs [EUR]: ", DSO_costs)
println("S_sub [kVA]:     ", S_sub)
println(" ")

fig_dict = ("p_self_sufficiency" => hist1,
    # "q_self_sufficiency" => hist2,
    "p_self_consumption" => hist3,
    # "q_self_consumption" => hist4,
    "user_costs" => hist5,
    "s_grid_max" => hist6,
    "p_imp_user" => hist7,
    "line_load" => hist8,
);
[savefig(value, file_name *"__"* key) for (key,value) in fig_dict]












# ###
# using Luxor
# lines = [file_name,
#     "DSO_costs [EUR]: $DSO_costs",
#     "S_sub [kVA]:     $S_sub"]
# # @png begin
# Drawing(200,200, :png)
# Luxor.textbox(lines, leading=15, Point(100,0),)
# finish()
# preview()

# Drawing(1000, 1000, :png)
# background("black")
# sethue("red")
# fontsize(50)
# Luxor.text("hello world")
# finish()
# preview()