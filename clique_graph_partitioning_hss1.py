from dwave.system import LeapHybridSampler
import dimod
import sys

# Graph partitioning on clique
N = int(sys.argv[1])
gamma = 3
linear = (N - 1) * (1 - gamma)
quad = (2 * gamma) - 2


Q = {}
for i in range(N):
    Q[i, i] = linear
    for j in range(i + 1, N):
        Q[i, j] = quad

chainstrength = 20

offset = gamma * N**2 / 4
bqm = dimod.BinaryQuadraticModel.from_qubo(Q, offset=offset)

sampler = LeapHybridSampler(profile='hss')
response = sampler.sample(bqm)

print(response.info)
for sample, energy in response.data(['sample', 'energy']):
    print(sample, energy)
