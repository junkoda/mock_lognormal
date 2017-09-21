# mock_lognormal
Lognormal mock generator (large-scale structure cosmology)

Example: 
mock_lognormal --np=100 --boxsize=1000 --z=1 --omega-m=0.308 planck_matterpower.dat

Input file:
planck_matterpower.dat: Tabulated k P(k) at z=0

--boxsize length of the preiodic box on a side
--nc      number of grids per dimension
--z       Redshift of the mock
--omega_m Omega matter at present
--np      Number of output particles

--lognormal outputs become lognormal density mock (default is Gaussian)

--rsd
Power spectrum becomes Kaiser
Ps(k, mu) = (1 + f*mu^2)^2 P(k)
where mu = kz/k

Default is real space without redshift space.

The particle distribution is random inside a cell of size boxsize/nc.

--fix-amplitude
Fix the amplitude of delta(k), and assign randomness only in the phase.

## Output

Cartesian corrdinate x, y, z in comoving coodinate.

Column 1: x [1/h Mpc]
Column 2: y [1/h Mpc]
Column 3: z [1/h Mpc]

Line of sight is z direction when --rsd is on.