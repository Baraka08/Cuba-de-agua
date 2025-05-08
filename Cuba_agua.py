# -*- coding: utf-8 -*-
"""
HOla que tal
241Am-Be point neutron source surrounded by air.

Emilio Castro González <emilio.castro@upm.es>
Nuria García Herranz <nuria.garcia.herranz@upm.es>

Geometry:
    Source is located in the center of a 10 m radius sphere filled with low density air or vacuum

Material:
    Highly diluted air with the following composition in weight fraction:
        N-14 - 0.755268
        O-16 - 0.231881
        Ar   - 0.012727
        C    - 0.000124
        Air density: 1.205E-50 g/cm3
        
"""

import math
import openmc

detector_distance = 100 # cm

########################
# Materials definition #
########################

all_materials = openmc.Materials()

mat_air = openmc.Material(name='Air')
mat_air.add_nuclide('N14', 0.755268, 'wo')
mat_air.add_nuclide('O16', 0.231881, 'wo')
mat_air.add_element('Ar', 0.012727, 'wo')
mat_air.add_element('C', 0.000124, 'wo')
mat_air.set_density('g/cm3', 1.205E-3)
mat_air.temperature = 293

all_materials.append(mat_air)


mat_lead = openmc.Material(name='Lead')
mat_lead.add_nuclide('Pb206', 0.241, 'ao')
mat_lead.add_nuclide('Pb207', 0.221, 'ao')
mat_lead.add_nuclide('Pb208', 0.538, 'ao')
mat_lead.set_density('g/cm3', 11.34)
mat_lead.temperature = 293

all_materials.append(mat_lead)

mat_concrete = openmc.Material(name='Concrete')
mat_concrete.add_nuclide('H1', 0.013407, 'ao')
mat_concrete.add_nuclide('C12', 0.0011030, 'ao')
mat_concrete.add_nuclide('O16', 0.043887, 'ao')
mat_concrete.add_nuclide('Al27', 0.0017971, 'ao')
mat_concrete.add_nuclide('Si28', 0.016123, 'ao')
mat_concrete.add_nuclide('Ca40', 0.001894, 'ao')
mat_concrete.add_nuclide('Fe56', 0.00033448, 'ao')
mat_concrete.set_density('g/cm3', 2.2)
mat_concrete.temperature = 293

all_materials.append(mat_concrete)

mat_Lucite = openmc.Material(name='Lucite')
mat_Lucite.add_element('C', 0.6, 'ao')
mat_Lucite.add_nuclide('H1', 0.08, 'ao')
mat_Lucite.add_nuclide('O16', 0.32, 'ao')
mat_Lucite.set_density('g/cm3', 0.92)
mat_Lucite.temperature = 293

all_materials.append(mat_Lucite)

mat_Cadmium_foil = openmc.Material(name='Cadmium_foil')
mat_Cadmium_foil.add_nuclide('Cd113', 1, 'ao')
mat_Cadmium_foil.set_density('g/cm3', 8.6)
mat_Cadmium_foil.temperature = 293

all_materials.append(mat_Cadmium_foil)

mat_Indium_foil = openmc.Material(name='Indium_foil')
mat_Indium_foil.add_nuclide('In115', 1, 'ao')
mat_Indium_foil.set_density('g/cm3', 7.22)
mat_Indium_foil.temperature = 293

all_materials.append(mat_Indium_foil)

all_materials.export_to_xml()

#######################
# Surface definitions #
#######################

surf_sphe_chamber_ext = openmc.Sphere(r=1000, boundary_type='vacuum')
surf_sphe_chamber_int = openmc.Sphere(r=50)


####################
# Cells definition #
####################

cell_chamber = openmc.Cell(name='Sphere surrounding the source')
cell_chamber.region = -surf_sphe_chamber_ext & +surf_sphe_chamber_int
cell_chamber.fill = mat_air

cell_shield = openmc.Cell(name='Shielding')
cell_shield.region = -surf_sphe_chamber_int
cell_shield.fill = mat_polyethilene

#############################
# Final geometry definition #
#############################

root_universe = openmc.Universe(cells=[cell_chamber, cell_shield])
geom = openmc.Geometry()
geom.root_universe = root_universe
geom.export_to_xml()

###########################
# Fixed source definition #
###########################

src = openmc.IndependentSource()
src.space = openmc.stats.Point((0,0,0))
# https://openmc.discourse.group/t/source-energy-distribution-for-neutron-of-different-energy-group/477
# https://openmc.discourse.group/t/modelling-am241-be-source/4710/1
# https://doi.org/10.1016/0020-708X(70)90066-9
# Las probabilidades deben ir divididas por el ancho del intervalo. En este caso, igual porque todos son iguales
import numpy as np
lstEnergy = np.linspace(0.2e6, 11.2e6, 56).tolist()
lstProbability = [0.028, 0.028, 0.024, 0.0205, 0.024, 0.028, 0.0168, 0.0182, 0.0178, 
                  0.0183, 0.0202, 0.0202, 0.0201, 0.0225, 0.0286, 0.0351, 0.0362, 0.0324, 
                  0.0296, 0.0284, 0.0277, 0.0283, 0.0301, 0.0286, 0.0311, 0.0295, 0.0265, 
                  0.0241, 0.0216, 0.0184, 0.0168, 0.0169, 0.0162, 0.0146, 0.0134, 0.0143, 
                  0.0159, 0.0166, 0.0171, 0.0162, 0.0134, 0.0102, 0.0073, 0.0048, 0.0036, 
                  0.0040, 0.0053, 0.0064, 0.0064, 0.0058, 0.0048, 0.0035, 0.0022, 0.0011,
                  0.0003, 0.0001]

src.energy = openmc.stats.Tabular(lstEnergy, lstProbability, 'histogram')
src.strength = 5.2E6 #n/s

########################
# Calculation settings #
########################

settings = openmc.Settings()
settings.source = src
settings.batches = 20
settings.particles = 1000000
settings.run_mode = 'fixed source'
settings.temperature = {'method': 'interpolation'}
settings.export_to_xml()

###########
# Tallies #
###########

tallies = openmc.Tallies()

tally_name = 'Dose rate'
dose_tally = openmc.Tally(name=tally_name)
mesh = openmc.SphericalMesh(r_grid=[100 - 0.5, 100 + 0.5])
mesh_filter = openmc.MeshFilter(mesh)
# Those lines are used to multiply flux by dose conversion factors from ICRP116
# https://openmc.discourse.group/t/openmc-dose-calculation-how-to-determine-optimal-settings-for-accuracy/4957/5
# https://docs.openmc.org/en/latest/pythonapi/generated/openmc.data.dose_coefficients.html
energy_bins_n, dose_coeffs_n = openmc.data.dose_coefficients(
    particle='neutron')
energy_function_filter = openmc.EnergyFunctionFilter(energy_bins_n,
                                                     dose_coeffs_n,
                                                     interpolation='log-log')
dose_tally.filters = [energy_function_filter, mesh_filter]
dose_tally.scores = ['flux']
tallies.append(dose_tally)

tallies.export_to_xml()

##############
# Run OpenMC #
##############

openmc.run()

##################
# Postprocessing #
##################

sp = openmc.StatePoint(f'statepoint.{settings.batches}.h5')

dose_mean = sp.get_tally(name=tally_name).mean.item()
dose_std_dev = sp.get_tally(name=tally_name).std_dev.item()
dose_rel_dev = dose_std_dev / dose_mean * 100.
vol_sphere = 4 / 3 * math.pi * (mesh.r_grid[1]**3 - mesh.r_grid[0]**3)
dose = dose_mean / vol_sphere 
# Convert from pSv/s to uSv/h
dose = dose / 1.0E6 * 3600

dose = round(dose, 2)
dose_rel_dev = round(dose_rel_dev, 2)

sp.close()

print(f'Dose rate is {dose} microSv/h with a standard deviation of {dose_rel_dev}%.')

