import argparse
from subprocess import run
import pandas as pd
from datetime import datetime

t_ejec = datetime.now().strftime("%m%d_%H%M")

def qdah(start, end, runs):
    for i in range(start, end):
        N = 4 * i**3
        B = [2**b for b in range(5, min(N,1024)) if 2**b<(min(N,1024)+1)]
        for bloque in B:
            run(['make', 'clean'])
            run(['make', f'CPPFLAGS=-DSIZE={N} -DBLOCK={bloque}'])
            for j in range(runs):
                run(['./Qdah'])
                readed = pd.read_csv("metrica.temp", header = 0)
                readed.insert(0, "BLOCK_SIZE", bloque)
                readed = readed.rename(index = {0: N})
                readed.to_csv(f"qdah_{t_ejec}.csv", index = True, header = False,
                                mode = 'a')
    run(['make', 'clean'])

def main():
    if args.qdah:
        print(args.qdah[0], args.qdah[1], args.qdah[2])
        d = [("N", "BLOCK_SIZE", "GFLOPS")]
        header = pd.DataFrame(data = d)
        header.to_csv(f"qdah_{t_ejec}.csv", index = False, header = False,
                        mode = 'w')
        qdah(args.qdah[0], args.qdah[1], args.qdah[2])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--qdah', nargs = 3, type = int,
                            help = "Corre el qdah. Pide 3 enteros: el m inicial \
                            (start), el m final más 1 (end) y la cantidad de \
                            veces que debe correr con cada m. La cantidad de \
                            partículas N viene dada por 4*m^3.")
    args = parser.parse_args()

    main()
