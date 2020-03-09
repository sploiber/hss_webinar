from minorminer import find_embedding
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import FixedEmbeddingComposite
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

dwave_sampler = DWaveSampler()
embedding = find_embedding(Q, dwave_sampler.edgelist)
chain_lengths = [len(embedding[k]) for k in embedding]
print("Total qubits ", sum(chain_lengths))
print("Longest chain ", max(chain_lengths))

sampler = FixedEmbeddingComposite(DWaveSampler(), embedding)
response = sampler.sample(bqm, chain_strength=chainstrength, num_reads=1000)

for sample, energy, num, cbf in response.data(['sample', 'energy', 'num_occurrences', 'chain_break_fraction']):
    print(sample, energy, num, cbf)
