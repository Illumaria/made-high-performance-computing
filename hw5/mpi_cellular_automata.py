#!/usr/bin/env python
import random
from argparse import ArgumentParser
from typing import List, NoReturn

from mpi4py import MPI

NUM_BORDERS = 2
LEFT_TO_RIGHT = 0
RIGHT_TO_LEFT = 1


def print_array(array: List[int]) -> NoReturn:
    print("".join(["⬜" if x else "⬛" for x in array]), flush=True)


def exchange_ghost_cells(
    array: List[int],
    comm: MPI.Intracomm,
    rank: int,
    world_size: int,
    periodic: bool = True,
) -> List[int]:
    if world_size > 1:
        # send leftmost real cell to the left neighbour rightmost ghost cell
        left_ghost_value = array[1] if periodic or 0 < rank < world_size else 0
        comm.send(left_ghost_value, dest=(rank - 1) % world_size, tag=RIGHT_TO_LEFT)
        # send rightmost real cell to the right neighbour leftmost ghost cell
        right_ghost_value = array[-2] if periodic or -1 < rank < world_size - 1 else 0
        comm.send(right_ghost_value, dest=(rank + 1) % world_size, tag=LEFT_TO_RIGHT)
        # receive rightmost ghost cell from the right neighbour
        array[-1] = comm.recv(tag=RIGHT_TO_LEFT)
        # receive leftmost ghost cell from the left neighbour
        array[0] = comm.recv(tag=LEFT_TO_RIGHT)
    else:
        array[-1] = array[1] if periodic else 0
        array[0] = array[-2] if periodic else 0

    return array


def generate_rule_patterns(rule_id: int) -> List[List[int]]:
    """
    For a given rule ID generates a list of patterns
    that correspond to 1's as new cell state.

    Example:
        >>> generate_rule_patterns(rule_id=110)
        [[0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]]
    """
    binary_rule_id = f"{rule_id:08b}"
    patterns = [f"{i:03b}" for i, c in enumerate(binary_rule_id[::-1]) if int(c) == 1]
    patterns = [[int(c) for c in pattern] for pattern in patterns]
    return patterns


def apply_rule(array: List[int], rule_patterns: List[List[int]]) -> List[int]:
    new_array = array.copy()
    for i in range(1, len(array) - 1):
        new_array[i] = int(array[i - 1:i + 2] in rule_patterns)
    return new_array


def run_cellular_automata(
    n_epochs: int, size: int, rule_id: int, periodic: bool, show_progress: bool
) -> NoReturn:
    start_time = MPI.Wtime()
    comm = MPI.COMM_WORLD
    world_size = comm.Get_size()
    rank = comm.Get_rank()

    if size % world_size != 0:
        MPI.Finalize()
        raise NotImplementedError(
            f"The program is implemented only for sizes that are multiples "
            f"of world size which is equal to {world_size}. "
            f"Please modify the 'size' accordingly (currently: {size})."
        )

    rank_arr_size = size // world_size
    array = [random.choice([0, 1]) for _ in range(rank_arr_size + NUM_BORDERS)]
    rule_patterns = generate_rule_patterns(rule_id)

    for _ in range(n_epochs):
        array = exchange_ghost_cells(array, comm, rank, world_size, periodic)
        full_array = comm.gather(array, root=0) if world_size > 1 else [array]
        if show_progress and rank == 0:
            print_array([elem for subarray in full_array for elem in subarray[1:-1]])

        # apply rule
        array = apply_rule(array, rule_patterns)

    full_array = comm.gather(array, root=0) if world_size > 1 else [array]
    if rank == 0:
        if show_progress:
            print_array([elem for subarray in full_array for elem in subarray[1:-1]])
        print(
            f"Time elapsed for world size of {world_size}: "
            f"{MPI.Wtime() - start_time} s"
        )

    MPI.Finalize()


def setup_parser(parser: ArgumentParser) -> NoReturn:
    parser.add_argument("--n-epochs", type=int, default=10)
    parser.add_argument("--size", type=int, default=10)
    parser.add_argument("--rule-id", type=int, default=110)
    parser.add_argument("--periodic", action="store_true")
    parser.add_argument("--show-progress", action="store_true")


def main() -> NoReturn:
    parser = ArgumentParser(
        prog="cellular_automata",
        description="MPI implementation of cellular automata",
    )
    setup_parser(parser)

    arguments = parser.parse_args()

    run_cellular_automata(
        arguments.n_epochs,
        arguments.size,
        arguments.rule_id,
        arguments.periodic,
        arguments.show_progress,
    )


if __name__ == "__main__":
    main()
