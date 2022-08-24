import scipy.constants as phyconst
import numpy as np

𝜀0 = phyconst.epsilon_0
me = phyconst.m_e
c = phyconst.c
qe = phyconst.e
pi = phyconst.pi
Lambda = 0.8e-6
Omega = 2*pi*c/Lambda

Main(
    geometry = '2Dcartesian',
    grid_length = [64*2*pi, 32*2*pi],
    cell_length = [0.1*2*pi,0.1*2*pi],
    simulation_time = 3*64*2*pi,
    timestep = 0.01*2*pi,
    number_of_patches = [32, 16],
    EM_boundary_conditions = [["silver-muller", "silver-muller"], ["silver-muller", "silver-muller"]],
    reference_angular_frequency_SI = Omega,
)
MovingWindow(
    time_start = 70*2*pi,
    velocity_x = 0.85,
)

def ND(x, y):
    if x>25*2*pi:
        if abs(y - 16*2*pi)<=15*2*pi:
            return 20
    elif x>=20 and x<=25:
        if abs(y - 16*2*pi)<=15*2*pi:
            return (x-20*2*pi)*4
    else:
        return 0

Species(
    name = 'electron',
    position_initialization = 'random',
    momentum_initialization = 'cold',
    particles_per_cell = 30,
    mass = 1,
    number_density = ND,
    charge = -1,
    mean_velocity = 0,
    temperature = 0,
    boundary_conditions = [['remove', 'remove'],['remove', 'remove']],
    radiation_model = "Monte-Carlo",
    radiation_photon_species = 'photon',
    radiation_photon_gamma_threshold = 1,
)
Species(
    name = 'proton',
    position_initialization = 'random',
    momentum_initialization = 'cold',
    particles_per_cell = 15,
    mass = 1836,
    number_density = ND,
    charge = 1,
    mean_velocity = 0,
    temperature = 0,
    boundary_conditions = [['remove', 'remove'],['remove', 'remove']],
)
Species(
    name = 'photon',
    position_initialization = 'random',
    momentum_initialization = 'cold',
    particles_per_cell = 0,
    mass = 0,
    number_density = 0,
    charge = 0,
    mean_velocity = 0,
    temperature = 0,
    boundary_conditions = [['remove', 'remove'],['remove', 'remove']],
)

LaserGaussian2D(
    box_side         = "xmin",
    a0               = 500,
    omega            = 1,
    focus            = [10*2*pi, 16*2*pi],
    waist            = 5*2*pi,
    incidence_angle  = 0.,
    polarization_phi = 0.,
    ellipticity      = 0.,
    time_envelope    = tgaussian(start=1*2*pi, duration=18*2*pi, center=10, order=2),
)

RadiationReaction(
  minimum_chi_continuous = 1e-3,
  minimum_chi_discontinuous = 1e-2,
  Niel_computation_method = "table",
)

DiagScalar(
    every = 100,
    vars = ["Utot", "Ukin", "Uelm", "ExMin", "ExMax"],
    precision = 10,
)
# 场分布
DiagFields(
    every = 200,
    time_average = 2,
    fields = ["Ex", "Ey", "Ez"],
)
DiagFields(
    every = 200,
    time_average = 2,
    fields = ["Bx", "By", "Bz"],
)
DiagFields(
    every = 200,
    time_average = 2,
    fields = ["Jx", "Jy", "Jz"],
)
# 密度分布
DiagParticleBinning(
    name = "electron density",
    deposited_quantity = "weight",
    every = 200,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["moving_x", 0., 64*2*pi, 640],
        ["y", 0, 32*2*pi, 320]
    ]
)
DiagParticleBinning(
    name = "proton density",
    deposited_quantity = "weight",
    every = 200,
    time_average = 1,
    species = ["proton"],
    axes = [
        ["moving_x", 0., 64*2*pi, 640],
        ["y", 0, 32*2*pi, 320]
    ]
)
DiagParticleBinning(
    name = "photon density",
    deposited_quantity = "weight",
    every = 200,
    time_average = 1,
    species = ["photon"],
    axes = [
        ["moving_x", 0., 64*2*pi, 640],
        ["y", 0, 32*2*pi, 320]
    ]
)
# 能谱分布
DiagParticleBinning(
    name = "log-electron",
    deposited_quantity = "weight",
    every = 200,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["gamma", 1, 5e4, 1000, "logscale"]
    ]
)
DiagParticleBinning(
    name = "log-photon",
    deposited_quantity = "weight",
    every = 200,
    time_average = 1,
    species = ["photon"],
    axes = [
        ["gamma", 1, 1e4, 1000, "logscale"]
    ]
)
DiagParticleBinning(
    name = "linear-electron",
    deposited_quantity = "weight",
    every = 200,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["gamma", 1, 5e4, 2000]
    ]
)
DiagParticleBinning(
    name = "linear-photon",
    deposited_quantity = "weight",
    every = 200,
    time_average = 1,
    species = ["photon"],
    axes = [
        ["gamma", 1, 1e4, 2000]
    ]
)
# 
DiagTrackParticles(
    species = "electron",
    every = 1000,
    attributes = ["x", "y", "px", "py", "chi"]
)