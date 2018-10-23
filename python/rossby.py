"""Nonlinear Rossby soliton test

Reference: Boyd (1980)

"""

import os
import sys
import time

import numpy as np
import dedalus.public as de
from dedalus.extras import flow_tools
from dedalus.tools import post

import logging
logger = logging.getLogger(__name__)


# Grid parameters
Lx = 40
Ly = 10

nx = 400
ny = 100

# Model parameters
A = 0.12
B = 0.394

x = de.Fourier('x',nx,  interval=[-Lx/2, Lx/2])
y = de.Chebyshev('y',ny, interval=[-Ly/2, Ly/2])

domain = de.Domain([x,y], grid_dtype='float')

problem = de.IVP(domain, variables=['u', 'v', 'phi'])

problem.parameters['Lx'] = Lx
problem.parameters['Ly'] = Ly
problem.parameters['a'] = 0
problem.parameters['β'] = 1
problem.substitutions['f'] = "β*y"
problem.substitutions['vol_avg(A)']   = 'integ(A)/(Lx*Ly)'

problem.add_equation("dt(u) + a*u + dx(phi) - f*v = -u*dx(u) - v*dy(u)")
problem.add_equation("dt(v) + a*v + dy(phi) + f*u = -u*dx(v) - v*dy(v)")
problem.add_equation("dt(phi) + dx(u) + dy(v) = -dx(u*phi) - dy(v*phi)")

problem.add_bc("left(v) = 0")
problem.add_bc("right(v) = 0")

# Build solver
ts = de.timesteppers.RK443
IVP = problem.build_solver(ts)

IVP.stop_sim_time = 10.
IVP.stop_wall_time = np.inf
IVP.stop_iteration = 10000000
dt = 0.025
logger.info('Solver built')

# Initial Conditions
eta = domain.new_field()
v0 = IVP.state['v']
u0 = IVP.state['u']
phi0 = IVP.state['phi']

xx, yy = domain.grids()

eta['g'] = A/np.cosh(B*xx)**2
deta = eta.differentiate('x')
u0['g'] = eta['g'] * (-9 + 6*yy**2)/4. * np.exp(-yy**2/2)
v0['g'] = 2 * deta['g'] * yy * np.exp(-yy**2/2)
phi0['g'] = eta['g'] * (3 + 6*yy**2)/4. * np.exp(-yy**2/2)

data_dir = 'rossby_soliton'
if domain.distributor.rank == 0:
    if not os.path.exists('{:s}/'.format(data_dir)):
        os.mkdir('{:s}/'.format(data_dir))

# Analysis
snapshots = IVP.evaluator.add_file_handler(os.path.join(data_dir,'snapshots'), sim_dt=1., max_writes=50)
snapshots.add_system(IVP.state)

timeseries = IVP.evaluator.add_file_handler(os.path.join(data_dir,'timeseries'),iter=10, max_writes=np.inf)
timeseries.add_task("vol_avg(sqrt(u*u))",name='urms')
timeseries.add_task("vol_avg(sqrt(v*v))",name='vrms')

analysis_tasks = [snapshots, timeseries]

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while IVP.ok:
        dt = IVP.step(dt)
        if (IVP.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(IVP.iteration, IVP.sim_time, dt))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %IVP.iteration)
    logger.info('Sim end time: %f' %IVP .sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

    for task in analysis_tasks:
        logger.info(task.base_path)
        post.merge_analysis(task.base_path)

