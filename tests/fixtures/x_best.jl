println("# Obtain the solution")
println("Retreiving solution...")
params = [kPhEnsa, kDpEnsa, kPhGw, kDpGw1, kDpGw2, kWee1, kWee2, kPhWee, kDpWee, kCdc25_1, kCdc25_2, kPhCdc25, kDpCdc25, kDipEB55, kAspEB55]
paramvalues = Dict()
for param in params
    paramvalues[param] = JuMP.value.(param)
end

paramvalues


