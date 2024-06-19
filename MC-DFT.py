#!/usr/bin/env python
import random, os, re, shutil, subprocess, math
import numpy as np
## routine to get final F (free energy including kinetic and electronic entopy)
def getF(file):
    output = subprocess.getoutput("grep 'sigma'  OUTCAR | tail -2")
    energy=[]
    for line in output.split('\n'):
        try:
                float(line.split()[-1])
                energy.append(line.split()[-1])
        except ValueError:
                print("no energy is obtained")
    if (float(energy[0])-float(energy[1])) > 1E-7:
       exit('Formation energy is not converged')
    return float(energy[1])

## routine to scale velocities by power of mass ratio
#def vScale(vel,ratio):
# power = 0, 0.5, 1 for constant velocity, kinetic energy or momentum
#    power = 0.5
#    vtmp = [float(vel[i])*ratio**power for i in range(3)]
#    return " ".join(format(v, ".6e") for v in vtmp)+'\n'


## preliminary setup
random.seed()
kB = 1./11604.

## Set Monte Carlo Temperature (K)
Temp = 100


## run initial step
subprocess.run([ 'mpirun', '-n','128', 'vasp_std'])

with open('INCAR', 'r') as f:
    incar = f.readlines(); f.close()

## initial energies used as first comparison 
Fi = getF("./OUTCAR")



for step in range(0, 1000):
        ## line of atom types
        with open('POSCAR', 'r') as f:
            poscar = f.readlines(); f.close()
        
        
        ## line of atom numbers each_num_atom =int(poscar[2].split()[0])/ntypes
        types = poscar[5].split()
        ntypes = len(types)
        atom_num = list(map(int,poscar[6].split()))
        natoms = sum(atom_num)
        start = [0]
        for t in range(1,ntypes):
            start.append(start[t-1]+atom_num[t-1])

        # get masses for velocity rescale
        #mass = []
        #with open('./POTCAR', 'r') as f:
        #    potcar = f.readlines(); f.close()
        #    for line in potcar:
        #        if "POMASS" in line: mass.append (float((line.split()[2]).replace(';', '')))

        # prepare subdirectory
        if not os.path.exists("step_"+str(step)):
           os.mkdir("step_"+str(step))
           shutil.copy("KPOINTS","step_"+str(step))
           shutil.copy("POTCAR","step_"+str(step))
           shutil.copy("INCAR","step_"+str(step))
        #for line in incar:
        #    if "NSW" in line: incar[incar.index(line)] = "NSW = 1\n"
        #    if "TEBEG" in line: Temp = float(line.split()[-1])

        #with open('./mcdir/INCAR', 'w') as f:
        #    f.writelines(incar); f.close()
        abspath = os.getcwd()
        dname = os.path.dirname(abspath)
        runpath = abspath+'/step_'+str(step)
        os.chdir(runpath)
        # create list of swaps
        swaps = []
        for t1 in range(0,ntypes):
            for t2 in range(t1+1,ntypes):
                swaps.append([t1, t2])

        nswaps = len(swaps)

        ## Create static run of swapped structure
        # choose swap
        swap = swaps[random.randint(0,nswaps-1)]
        swap0 = random.randint(0,atom_num[swap[0]]-1); swap1 = random.randint(0,atom_num[swap[1]]-1)
        atom0 = start[swap[0]]+swap0; atom1 = start[swap[1]]+swap1
        print("attempt t0=%d a0=%d with t1=%d a1=%d" % (swap[0],swap0,swap[1],swap1)),

        # POSCAR offset for first atom is 8 lines
        atom0 = atom0+8; atom1=atom1+8
        # swap coordinates
        poscar[atom0], poscar[atom1] = poscar[atom1], poscar[atom0]
        #coord_a = poscar[atom0].split()
        #coord_b = poscar[atom1].split()
        #poscar[atom0] = " %s  %s    %s %s %s 0   0  0 \n" % (coord_a[0],coord_a[1],coord_b[2],coord_b[3],coord_b[4])
        #poscar[atom1] = " %s  %s    %s %s %s 0   0  0 \n" % (coord_b[0],coord_b[1],coord_a[2],coord_a[3],coord_a[4])
        #poscar[atom1] = [coord_b[0],coord_b[1],coord_a[2],coord_a[3],coord_a[4],coord_a[5],coord_a[6],coord_a[7]]
        ## Swapping velocities is allowed by unnecessary
        ## Uncomment following lines to swap and scale velocities if desired
        #v0_new = vScale(poscar[atom1+natoms+1].split(),mass[swap[1]]/mass[swap[0]])
        #v1_new = vScale(poscar[atom0+natoms+1].split(),mass[swap[0]]/mass[swap[1]])
        #poscar[atom0] = v0_new
        #poscar[atom1] = v1_new
        #incar[7] = 'read_data    %s\n' % ('structure.lammps_vs'+str(step))
        #incar[-1] = 'write_data    %s\n' % ('structure.lammps_vs'+str(step+1))
        # omit predictors and correctors since they are incorrect
        #with open("./structure.lammps_MID"+str(step+1), 'w') as f:
        #    f.writelines(poscar[0:20+natoms]); f.close()
        with open("POSCAR", 'w') as f:
             f.writelines(poscar); f.close()


        #abspath = os.path.abspath(__file__)
        #abspath = os.getcwd()
        #dname = os.path.dirname(abspath)
        #runpath = dname+'/step_'+step
        #os.chdir(runpath)

        ## Run VASP
        subprocess.run([ 'mpirun', '-n','128', 'vasp_std'])
        #subprocess.run('lmp -in lammps.in>lammps.out', shell=True)


        ## find dF and accept or reject
        Ff = getF("OUTCAR")


        dF=Ff-Fi
        print("Ff=%.3f Fi=%.3f dF=%.3f" % (Ff, Fi, dF)),
        if dF < 0:
            prob = 1.000; r=1.000; accept = 1
        else:
            prob = math.exp(-dF/(kB*Temp))
            r = random.random()
            if prob < r:
                accept = 0
                print("prob=%.3f rand=%.3f reject" % (prob,r)) #os.remove("mcdir/WAVECAR")
            else:
                accept = 1
        with open("../energy.txt", "a") as text_file:
            text_file.write(str(step)+"    "+str(Ff)+"     "+ str(accept)+"\n")

        if accept == 1:
            Fi=Ff
            print("prob=%.3f rand=%.3f accept" % (prob,r))
            shutil.copy("CONTCAR","../POSCAR")

        
        os.chdir(abspath)
