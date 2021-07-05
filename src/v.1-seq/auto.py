import os
import argparse
from subprocess import run
import pandas as pd
import matplotlib.pyplot as plt



def test(comp):
    for i in comp:
        run(['make', 'clean'])
        run(['make','CC='+i[0] ,'CFLAGS='+i[1]])
        for j in range(10):
            run(['./tiny_md'])
            readed = pd.read_csv("metrica.temp",header=0)
            readed.insert(0,"Flags","".join(i[1:]))
            readed=readed.rename(index={0: i[0]})
            readed.to_csv("medidas.data",index=True,header=False,mode='a')


def plot():
    return 0

def main():    
    if args.add:
        arg = args.add
        d = {'comp': [arg[0]], 'opts': [" ".join(arg[1:])]}
        new = pd.DataFrame(data=d)
        new.to_csv("barrido.txt",index=False,header=False,mode='a')
    if args.compile:
        d = [("Flags","Gflops_Force")]
        header = pd.DataFrame(data=d)
        header.to_csv("medidas.data",index=False,header=False,mode='w')
        comps = pd.read_csv("barrido.txt",header=0)
        test(comps.values)
    if args.plot:
        data = pd.read_csv("medidas.data",header=0,index_col=0)
        index = data.index.drop_duplicates()
        for i in index:
            promedios = pd.DataFrame(columns=data.columns[1:])
            for j in data.loc[i]['Flags'].drop_duplicates():
                promedios.loc[j] = data.loc[i][data.loc[i]['Flags'] == j].mean()
            promedios.plot.bar()
            plt.suptitle(i)
            plt.savefig('plot'+i+'.png', dpi=300, bbox_inches='tight')
            #plt.show()
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--plot',help="Plotea un gráfico del rendimiento de las corridas",action="store_true")
    parser.add_argument('-c','--compile',help="Compila y corre el proyecto con distintas optimizaciones y/o compiladores pre definidos",action="store_true")
    parser.add_argument('-a','--add',help="Agrega la compilacion", nargs=2,metavar=('CC','"FLAGS"'))
    args = parser.parse_args()   

    main()

