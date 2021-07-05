import argparse
from subprocess import run
import pandas as pd
from datetime import datetime

t_ejec = datetime.now().strftime("%m%d_%H%M")

def compiladores(comp):
    for i in comp:
        run(['make', 'clean'])
        run(['make', 'CC='+i[0], 'CFLAGS='+i[1]])
        for j in range(10):
            run(['./tiny_md'])
            readed = pd.read_csv("metrica.temp", header = 0)
            readed.insert(0, "Flags", "".join(i[1:]))
            readed = readed.rename(index = {0: i[0]})
            readed.to_csv(f"comp_{t_ejec}.csv", index = True, header = False,
                            mode = 'a')
    run(['make', 'clean'])

def particulas(start, end, runs):
    for i in range(start, end):
        N = 4 * i**3
        run(['make', 'clean'])
        run(['make', f'CPPFLAGS=-DN={N}'])
        for j in range(runs):
            run(['./tiny_md'])
            readed = pd.read_csv("metrica.temp", header = 0)
            readed = readed.rename(index = {0: N})
            readed.to_csv(f"part_{t_ejec}.csv", index = True, header = False,
                            mode = 'a')
    run(['make', 'clean'])

def tiling(start, end, runs):
    for i in range(start, end):
        N = 4 * i**3
        B = [2**b for b in range(5, N) if 2**b<N+1]
        for bloque in B:
            run(['make', 'clean'])
            run(['make', f'CPPFLAGS=-DN={N} -DBLOCK={bloque}'])
            for j in range(runs):
                run(['./tiny_md'])
                readed = pd.read_csv("metrica.temp", header = 0)
                readed.insert(0, "BLOCK", bloque)
                readed = readed.rename(index = {0: N})
                readed.to_csv(f"tiling_{t_ejec}.csv", index = True, header = False,
                                mode = 'a')
    run(['make', 'clean'])

def main():
    if args.compiladores:
        d = [("Flags", "GFLOPS")]
        header = pd.DataFrame(data = d)
        header.to_csv(f"comp_{t_ejec}.csv", index = False, header = False,
                        mode = 'w')
        comps = pd.read_csv("comp.in", header=0)
        compiladores(comps.values)
    if args.particulas:
        print(args.particulas[0], args.particulas[1], args.particulas[2])
        d = [("N", "GFLOPS")]
        header = pd.DataFrame(data = d)
        header.to_csv(f"part_{t_ejec}.csv", index = False, header = False,
                        mode = 'w')
        particulas(args.particulas[0], args.particulas[1], args.particulas[2])
    if args.tiling:
        print(args.tiling[0], args.tiling[1], args.tiling[2])
        d = [("N", "BLOCK", "GFLOPS")]
        header = pd.DataFrame(data = d)
        header.to_csv(f"tiling_{t_ejec}.csv", index = False, header = False,
                        mode = 'w')
        tiling(args.tiling[0], args.tiling[1], args.tiling[2])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--compiladores', action = "store_true",
                            help = "Compila y corre el proyecto con distintas \
                            optimizaciones y/o compiladores pre definidos en \
                            el archivo 'comp.in'.")
    parser.add_argument('-p', '--particulas', nargs = 3, type = int,
                            help = "Compila y corre el proyecto con un dado \
                            compilador y ciertas optimizaciones, barriendo la \
                            cantidad de partículas. Se deben dar 3 enteros: el \
                            m inicial (start), el m final más 1 (end) y la \
                            cantidad de veces que debe correr con cada m. \
                            La cantidad de partículas N viene dada por 4*m^3.")
    parser.add_argument('-t', '--tiling', nargs = 3, type = int,
                            help = "Compila y corre el proyecto con un dado \
                            compilador y ciertas optimizaciones, barriendo la \
                            cantidad de partículas y, por cada una de ellas, \
                            va cambiando el tamaño de los bloques según 2^i. \
                            Se deben dar 3 enteros: el m inicial (start), el m \
                            final más 1 (end) y la cantidad de veces que debe \
                            correr con cada m. La cantidad de partículas N \
                            viene dada por 4*m^3.")
    args = parser.parse_args()

    main()
