
import numpy as np
import os
import sys


# main(R)
# pots = ['D2ClNaPerturbTrain7r11','PBED2NaClTrain1AllTraining']
# pot = 'D2ClNaPerturbTrain7r11'

if len(sys.argv) < 2:
    print('Usage: mdext...py <pot file stem wo ext.>')
    exit(1)

pot = str(sys.argv[1])

# os.mkdir(pot)
# os.chdir(pot)
# run cavitation sims
# batch them 
#Rs = np.arange(0.2,9.1,0.1) 
#Rs = np.arange(0.0,0.3,0.1)

# Rs = np.arange(0.0,9.1,0.1)
# LO 5.0
# Rs = np.arange(5.1,9.1,0.1)
LO = 6.9  # what to run next
Rs = np.arange(LO,9.1,0.1)
Rs = list(Rs)
Rs = [round(R,1) for R in Rs]
#print(Rs)
#exit()
for i, R in enumerate(Rs):
    #print(R)
    #break
    R = round(R,1)
    print(f"launching cavity {R}")
    LOlast = round(LO-0.1,1)
    # if R == 0.0:  # initial run
    #     startfile = 'liquid.data'
    # if R == 5.1:  # continue initial run
    #    startfile = '5.0.cavity.data'
    if R == LO:  # continue initial run
       startfile = f'{LOlast}.cavity.data'
    else:  # normal sequence within a run
        prev_R = Rs[Rs.index(R)-1]
        startfile = str(prev_R)+'.cavity.data'  
    randomSeed = np.random.randint(0,1000)
    # switch if dpmd or cavity
    # pot='PBED2NaClTrain1AllTraining'
    # os.system("nohup bash cavityDPMD.job %s %s %s %s > LOG &"%(str(R),startfile,str(randomSeed),pot))
    # os.system("nohup bash cavity.job %s %s %s > LOG &"%(str(R),startfile,str(randomSeed)))
    os.system(f'python ../mdextCavitationDPMD.py {R} {startfile} {randomSeed} {pot}')
    # break
    # currdatafile = str(R) + '.cavity.data'
    # pwd = os.getcwd()
    # while not (os.path.exists(pwd+'/'+currdatafile)):
        # time.sleep(1)   
    # if i==1: break  # HACK
print('Done!')