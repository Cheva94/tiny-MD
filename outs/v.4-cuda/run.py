import sys
import getopt
from json import load, dumps, JSONDecodeError
from target import Target
from subprocess import run
from os import path as ph, chdir as ch, environ


""" Run variables """
CORE = 1                # Selected core[s] to run the likwid-perfctr.
COLLECT = False         # Collect perf data.

""" Counter groups, and collect metadata """
FL_SP = ('FLOPS_SP', ['Runtime (RDTSC) [s]',
                      'SP [MFLOP/s]'])

BR = ('BRANCH', ['Branch rate',
                 'Branch misprediction rate',
                 'Branch misprediction ratio',
                 'Instructions per branch'])

L2 = ('L2CACHE', ['L2 request rate',
                  'L2 miss rate',
                  'L2 miss ratio'])

""" Defined targets """
targets = [
    Target(name='T_CLEAN', stats_collectors=[FL_SP], flags=['-O0']),
    Target(name='T_O1', stats_collectors=[FL_SP], flags=['-O1', '-march=native']),
    Target(name='T_O2', stats_collectors=[FL_SP], flags=['-O2', '-march=native']),
    Target(name='T_O3', stats_collectors=[FL_SP], flags=['-O3', '-march=native']),
    Target(name='T_O2', stats_collectors=[FL_SP], flags=['-O3+', '-march=native']),
    Target(name='T_O2', stats_collectors=[FL_SP], flags=['-O3+', '-funroll-loops']),
    Target(name='T_Of', stats_collectors=[FL_SP], flags=['-Ofast'])
]

""" Directory handlers """
def fw(n): return ch(ph.join(ph.curdir, n))
def bw(): return ch('../')


def runner(t: Target, log_file):
    # Get env
    env = dict(environ)

    for ev in t.evars:
        env[ev[0]] = ev[1]

    fw(t.name)

    run_file = open('run.out', mode='w' if ph.isfile('run.out') else 'x')
    stats_file = open(
        'stats.json', mode='w' if ph.isfile('stats.json') else 'x')

    compile_res = run(['meson', 'compile'], shell=False,
                      capture_output=True, env=env)

    if(compile_res.returncode == 0):
        print(f'@ Compiled {t.name}\n@ Running ./headless')
        run_res = None
        if(t.stats_collectors is not None and COLLECT):
            for sc in t.stats_collectors:
                run_res = run(['likwid-perfctr', '-C', f'{CORE}', '-g', sc[0], '-m', '-O',
                              './headless'], shell=False, capture_output=True, env=env)
                if(run_res.returncode == 0):
                    output = run_res.stdout.decode('ascii')
                    run_file.write(output)
                    collect_res = collect(sc, output, stats_file, log_file)
                    if(collect_res is not None):
                        return collect_res
                else:
                    return f'@! Error running {t.name}/headless (likwid-perf)'

            bw()
            print(f'@ Collected performance stats: {t.name}/stats')
            return None

        else:
            run_res = run('./headless', shell=False,
                          capture_output=True, env=env)

            if(run_res.returncode == 0):
                bw()
                return None

            return f'@! Error running {t.name}/headless'

    else:
        bw()
        log_file.write(compile_res.stdout.decode('ascii'))
        return f'@! Error compiling {t.name}. See run.log for more details.'


def collect(sc, output, stats_file, log_file):
    print(f'@ Collecting data for group: {sc[0]}')

    try:
        lines = output.splitlines(keepends=False)

        regions_ids = []

        for id in range(0, len(lines)):
            line = lines[id]
            if(line.startswith('TABLE') and any(map(lambda x: 'Metric' in x, line.split(',')))):
                regions_ids.append(id)

        stats = {}
        for id in regions_ids:
            region_name = lines[id].split(',')[1].split(' ')[1]
            r_id = 1
            region_line = lines[id+r_id].split(',')
            stats[region_name] = {}
            while(region_line[0] != 'TABLE' and id+r_id < len(lines) - 1):
                if(region_line[0] in sc[1]):
                    region_stats = stats[region_name]
                    region_stats[region_line[0]] = region_line[1]
                r_id += 1
                region_line = lines[id+r_id].split(',')

        stats_file.write(dumps(stats))

        return None

    except Exception as err:
        log_file.write(str(err))
        return '@! Error collecting data. See run.log for more details.'


def configure(t, log_file):
    cmd = ['meson', 'setup', f'{t.name}']
    wipe_cmd = ['meson', 'configure']

    # Get env
    env = dict(environ)

    # Set compiler.
    if t.comp is not None:
        env['CC'] = t.comp

    # Set compiler flags.
    # c_args = list(map(lambda x: f'-Dc_args={x}', t.flags))
    c_args = [f"-Dc_args={','.join(t.flags)}"]
    print(c_args)

    if(t.stats_collectors is not None):
        c_args.append('-Dc_args=-DLIKWID_PERFMON')
        

    if(ph.isdir(ph.join(ph.curdir, t.name))):
        print('@ Target already exists, reconfiguring.')

        fw(t.name)
        
        res = run(wipe_cmd + c_args, shell=False, capture_output=True)

        log_file.write(res.stdout.decode('ascii'))

        bw()

        if(res.returncode != 0):
            return f'@! Error reconfiguring {t.name}. See run.log for more details.'

        print('@ Reconfigured succesfully.')

    else:
        print('@ Creating new target.')
        res = run(cmd + c_args, shell=False, capture_output=True, env=env)
        log_file.write(res.stdout.decode('ascii'))
        if(res.returncode != 0):
            return f'@! Error creating {t.name}. See run.log for more details.'

        print('@ Target created succesfully.')


def automatize(tgs):
    log_file = open('run.log', mode='w' if ph.isfile('run.log') else 'x')

    for tg in tgs:
        conf = configure(tg, log_file)
        if(conf is not None):
            print(conf)
            return
        runn = runner(tg, log_file)
        if(runn is not None):
            print(runn)
            return


def usage():
    """ Prints command usage. """
    print('Navier-Stokes automatization helper.')
    print('List of avaiable options:')
    print('-t TARGETNAME                  Automatize specific target.')
    print('-c CORE | CORE1 - COREN        Set cores to run performance counters.')
    print('-d                             Enable performance counters and collect data.')
    print('-s                             Generate summary from collected data.')


def clean():

    if(run('rm T_* -r', shell=True).returncode == 0):
        return '@ Build directories cleaned succesfully.'
    else:
        return '@ Error while cleaning directories.'


def summary(targets: list):

    summary_file = open('sum.json', mode='w' if ph.isfile('sum.json') else 'x')

    output = {}

    for t in targets:
        fw(t.name)
        s = open('stats.json', mode='r')
        try:
            sum = load(s)
            output[t.name] = sum
        except JSONDecodeError:
            output[t.name] = ''
        bw()
    summary_file.write(dumps(output))


def main():
    """ Main """
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "hsCdc:t:", ["help", "output="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    target = None
    clean_dirs = False
    s = False

    for o, a in opts:
        if o == '-c':
            global CORE
            CORE = a
        if o == '-t':
            if(a not in list(map(lambda x: x.name, targets))):
                print('Target not defined.')
                return
            target = [next((x for x in targets if x.name == a), None)]
        if o == '-d':
            global COLLECT
            COLLECT = True

        if o == '-C':
            clean_dirs = True

        if o == '-s':
            """ Summary """
            s = True

        elif o == '-h':
            usage()
            return

    if(s):
        summary(targets)
        return

    if(clean_dirs):
        print(clean())
        return

    automatize(targets if target is None else target)


if __name__ == '__main__':
    main()
