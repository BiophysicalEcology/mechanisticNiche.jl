
using CSV, DataFrames
using Plots

# define a function to plot some dataframes from micro
function plot_df(df::DataFrame)
    P = plot()
    for i in names(df)[3:end]
        plot!(df[:, i])
    end
    display(P)
end

metout = CSV.read("metout.csv", DataFrame, normalizenames=true)
soil = CSV.read("soil.csv", DataFrame, normalizenames=true)
soilmoist = CSV.read("soilmoist.csv", DataFrame, normalizenames=true)
soilpot = CSV.read("soilpot.csv", DataFrame, normalizenames=true)
tcond = CSV.read("tcond.csv", DataFrame, normalizenames=true)
sunsnow = CSV.read("sunsnow.csv", DataFrame, normalizenames=true)
# # print(metout)
# colNames = names(metout)
# colNames = names(soil)

# metout[:, [:JULDAY, :TIME]]

plot_df(soil)
plot!(metout[:, :TALOC], label="Taloc", linewidth=2)

plot_df(soilmoist)
plot_df(soilpot)
plot_df(tcond)
# plot_df(sunsnow)




