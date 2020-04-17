variables = [ :Cb
 :pCb
 :Wee
 :pWee
 :Cdc25
 :pCdc25
 :Gw
 :pGw
 :Ensa
 :pEnsa
 :pEB55
 :B55
]
variablevalues = Dict()
for v in variables
    variablevalues[string(v)] = Vector(JuMP.value.(eval(v)))
end

variablevalues


