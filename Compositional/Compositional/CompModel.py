import Flash
from Mixture import Mixture


psi = 3000
rankine = 160+460
molecules = ['C1', 'C4', 'C10']
zFracs = [0.6301, 0.1055, 0.2644]

mix = Mixture()
mix.initialize_fluids(molecules, zFracs)

Flash.flash_mixture(mix, psi, rankine)

print
print '-----Results-----'
print mix.zFactors
print mix.xFracs
print mix.yFracs
