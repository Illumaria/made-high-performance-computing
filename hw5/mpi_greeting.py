#!/usr/bin/env python
import random

from mpi4py import MPI


def main():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    available_ranks = list(range(size))
    processor_name = f"processor_{rank}"

    if rank == 0:
        greetings_chain = [(processor_name, rank)]
        dest = random.choice(available_ranks[1:])
        comm.ssend(greetings_chain, dest=dest)
    else:
        greetings_chain = comm.recv()
        print(f"Hello to {processor_name} from {greetings_chain[-1][0]}", flush=True)

        used_ranks = [x[1] for x in greetings_chain] + [rank]
        not_used_ranks = [x for x in available_ranks if x not in used_ranks]
        if len(not_used_ranks) > 0:
            greetings_chain.append((processor_name, rank))
            dest = random.choice(not_used_ranks)
            comm.ssend(greetings_chain, dest=dest)

    MPI.Finalize()


if __name__ == "__main__":
    main()
