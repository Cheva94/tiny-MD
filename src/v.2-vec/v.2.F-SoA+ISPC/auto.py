# import os
import argparse
from subprocess import run
import pandas as pd
# import matplotlib.pyplot as plt

def compiladores(comp):
    for i in comp:
        run(['make', 'clean'])
        run(['make', 'CC='+i[0], 'CFLAGS='+i[1]])
        for j in range(10):
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>> {i[0]}{i[1]}')
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Corrida {j+1} de 10')
            run(['./tiny_md'])
            readed = pd.read_csv("metrica.temp", header = 0)
            readed.insert(0, "Flags", "".join(i[1:]))
            readed = readed.rename(index = {0: i[0]})
            readed.to_csv("compiladores.data", index = True, header = False,
                            mode = 'a')
    run(['make', 'clean'])

def particulas(start, end):
    for i in range(start,end):
        N = 4 * (i+1)**3
        run(['make', 'clean'])
        run(['make', f'CPPFLAGS=-DN={N}'])
        for j in range(10):
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Corriendo con {N} particulas')
            print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Corrida {j+1} de 10')
            run(['./tiny_md'])
            readed = pd.read_csv("metrica.temp", header = 0)
            readed = readed.rename(index = {0: N})
            readed.to_csv("particulas.data", index = True, header = False,
                            mode = 'a')
    run(['make', 'clean'])

# def plot():
#     return 0

def main():
    if args.compiladores:
        d = [("Comp", "Flags", "GFLOPS")]
        header = pd.DataFrame(data = d)
        header.to_csv("compiladores.data", index = False, header = False,
                        mode = 'w')
        comps = pd.read_csv("compiladores.txt", header=0)
        compiladores(comps.values)
    if args.particulas:
        print(args.particulas[0],args.particulas[1])
        d = [("N", "GFLOPS")]
        header = pd.DataFrame(data = d)
        header.to_csv("particulas.data", index = False, header = False,
                        mode = 'w')
        particulas(args.particulas[0],args.particulas[1])
    # if args.add:
    #     arg = args.add
    #     d = {'comp': [arg[0]], 'opts': [" ".join(arg[1:])]}
    #     new = pd.DataFrame(data=d)
    #     new.to_csv("compiladores.txt",index=False,header=False,mode='a')
    # if args.plot:
    #     data = pd.read_csv("compiladores.data",header=0,index_col=0)
    #     index = data.index.drop_duplicates()
    #     for i in index:
    #         promedios = pd.DataFrame(columns=data.columns[1:])
    #         for j in data.loc[i]['Flags'].drop_duplicates():
    #             promedios.loc[j] = data.loc[i][data.loc[i]['Flags'] == j].mean()
    #         promedios.plot.bar()
    #         plt.suptitle(i)
    #         plt.savefig('plot'+i+'.png', dpi=300, bbox_inches='tight')
    #         plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--compiladores', action = "store_true",
                        help = "Compila y corre el proyecto con distintas \
                        optimizaciones y/o compiladores pre definidos")
    parser.add_argument('-p', '--particulas', nargs = 2, type=int, 
                        help = "Compila y corre el proyecto con un dado \
                        compilador y ciertas optimizaciones, barriendo la \
                        cantidad de partículas")
    # parser.add_argument('-p','--plot',help="Plotea un gráfico del rendimiento de las corridas",action="store_true")
    # parser.add_argument('-a','--add',help="Agrega la compilacion", nargs=2,metavar=('CC','"FLAGS"'))
    args = parser.parse_args()

    main()

