import os

os.system('cp ../06* .')

# Comment out old seeds already launched (to prevent relaunching), add more as needed
#seeds = range(0, 192)  # launched on cpu
#seeds = range(192, 384) # launched on cpu
#seeds = range(384, 640)  # launched on phi
seeds = [1]

for seed in seeds:
    print(f'Seed {seed}: ', end='', flush=True)
    runPath = f'Runs/seed{seed:04}'
    os.system(f'mkdir -p {runPath}')
    os.chdir(runPath)  
    os.system(f'sbatch ../../06run.job {seed}')
    os.chdir('../..')
