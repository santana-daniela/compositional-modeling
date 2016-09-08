from TernaryModel import create_ternary

psi = 3000 #Pressure in PSI
rankine = 160+460 #Temperature in Rankine
molecules = ['C1', 'C4', 'C10'] #Components
#specific_zFracs = [ 0.1055, 0.6301, 0.2644] # i and j are switched
amount = 10 #Number of test points in ternary plot
show_heat = True #Colorizes the plot based on fluid/liquid state.
show_tie = False #Create tie-lines
attempt_curve = False #Create boundary line around two-phase region
ternary_plot = create_ternary(molecules, psi, rankine, amt=amount, heatmap=show_heat, tie_lines=show_tie,
                              specific_zFracs=[], finish_curve=attempt_curve)