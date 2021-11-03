#!/usr/bin/env python
"""Python code to implement Conway's Game Of Life"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI

ON = 255
OFF = 0
NUM_BORDERS = 2
TOP_TO_BOTTOM = 0
BOTTOM_TO_TOP = 1
LEFT_TO_RIGHT = 2
RIGHT_TO_LEFT = 3
GOSPER_SIZE = 38


def random_grid(m: int, n: int, fraction_alive: float = 0.2) -> np.array:
    """Return a grid of [m x n] random values"""
    return np.random.choice([ON, OFF], m * n, p=[fraction_alive, 1 - fraction_alive]).reshape(m, n)


def gosper_glider_gun(m: int, n: int) -> np.array:
    """Return a grid in form of gosper glider gun"""
    glider_gun = np.array(
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
         [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
    ) * ON

    grid = np.zeros((m, n))
    start_row, start_col = (m - glider_gun.shape[0]) // 2, (n - glider_gun.shape[1]) // 2
    grid[start_row:start_row + glider_gun.shape[0], start_col:start_col + glider_gun.shape[1]] = glider_gun
    return grid


def apply_rules(grid: np.array) -> np.array:
    # copy grid since we require 8 neighbors
    # for calculation and we go line by line
    new_grid = grid.copy()
    m, n = grid.shape
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            # compute 8-neighbor sum
            total = int((grid[i, (j - 1) % n] + grid[i, (j + 1) % n] +
                         grid[(i - 1) % m, j] + grid[(i + 1) % m, j] +
                         grid[(i - 1) % m, (j - 1) % n] + grid[(i - 1) % m, (j + 1) % n] +
                         grid[(i + 1) % m, (j - 1) % n] + grid[(i + 1) % m, (j + 1) % n]) / ON)

            # apply Conway's rules
            if grid[i, j] == ON:
                if (total < 2) or (total > 3):
                    new_grid[i, j] = OFF
            else:
                if total == 3:
                    new_grid[i, j] = ON

    return new_grid


def exchange_ghost_cells(
        array: np.array,
        comm: MPI.Intracomm,
        rank: int,
        world_size: int
) -> np.array:
    if world_size > 1:
        # send the top row of real cells to the above neighbour bottom row of ghost cells
        comm.send(array[1, :], dest=(rank - 1) % world_size, tag=BOTTOM_TO_TOP)

        # send the bottom row of real cells to the below neighbour top row of ghost cells
        comm.send(array[-2, :], dest=(rank + 1) % world_size, tag=TOP_TO_BOTTOM)

        # send the leftmost column of real cells to the right column of ghost cells
        array[:, -1] = array[:, 1]

        # send the rightmost column of real cells to the left column of ghost cells
        array[:, 0] = array[:, -2]

        # receive rightmost ghost cell from the right neighbour
        array[-1, :] = comm.recv(tag=BOTTOM_TO_TOP)

        # receive leftmost ghost cell from the left neighbour
        array[0, :] = comm.recv(tag=TOP_TO_BOTTOM)
    else:
        array[-1, :] = array[1, :]
        array[0, :] = array[-2, :]
        array[:, -1] = array[:, 1]
        array[:, 0] = array[:, -2]

    return array


def run_game_of_life(n: int, fraction_alive: float, frame_rate: int, gosper: bool):
    comm = MPI.COMM_WORLD
    world_size = comm.Get_size()
    rank = comm.Get_rank()

    if world_size > n:
        raise NotImplementedError(f"The program is not implemented for world sizes "
                                  f"greater than the grid size. "
                                  f"Current grid size: {n}, world size: {world_size}.")

    rank_array_height = len(np.array_split(range(n), world_size)[rank]) + NUM_BORDERS
    rank_array_width = n + NUM_BORDERS
    if gosper and n >= GOSPER_SIZE:
        if rank == 0:
            rank_grid = gosper_glider_gun(rank_array_height, rank_array_width)
        else:
            rank_grid = np.zeros((rank_array_height, rank_array_width))
    else:
        if n < GOSPER_SIZE and rank == 0:
            print(f"Gosper glider gun requires grid of size {GOSPER_SIZE} or more. Using random grid instead.",
                  flush=True)
        rank_grid = random_grid(rank_array_height, rank_array_width, fraction_alive)

    full_grid = comm.gather(rank_grid, root=0) if world_size > 1 else [rank_grid]
    if rank == 0:
        full_grid = np.vstack([grid[1:-1] for grid in full_grid])
        fig, ax = plt.subplots()
        ax.set_xlim(left=0, right=full_grid.shape[0])
        img = ax.imshow(full_grid, interpolation="nearest")

    try:
        while True:
            rank_grid = exchange_ghost_cells(rank_grid, comm, rank, world_size)
            rank_grid = apply_rules(rank_grid)

            full_grid = comm.gather(rank_grid, root=0) if world_size > 1 else [rank_grid]
            if rank == 0:
                full_grid = np.vstack([grid[1:-1] for grid in full_grid])
                img.set_data(full_grid)
                fig.canvas.draw_idle()
                plt.pause(1 / frame_rate)
    except KeyboardInterrupt:
        MPI.Finalize()


def setup_parser(parser):
    parser.add_argument("--grid-size", type=int, default=100)
    parser.add_argument("--frame-rate", type=int, default=30)
    parser.add_argument("--fraction-alive", type=float, default=0.2)
    parser.add_argument("--gosper", action="store_true")


def main():
    parser = argparse.ArgumentParser(description="Runs Conway's Game of Life simulation")
    setup_parser(parser)
    args = parser.parse_args()

    run_game_of_life(args.grid_size, args.fraction_alive, args.frame_rate, args.gosper)


if __name__ == "__main__":
    main()
