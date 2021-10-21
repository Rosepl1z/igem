<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Document</title>
</head>
<body>
    <div class="pageName">
        <h1>gromacs.pl</h1>
        <div class="divider"></div>
    </div>
    <span>
        #! /usr/bin/perl -w

        print "Please input your PDB file(without'.pdb'):";
        
        $filename = <$ARGV[0]>;
        
        chomp $filename;
        
        print "$filename\n";
        
        system("gmx pdb2gmx -f $filename.pdb -o $filename\_processed.gro -water spce <span><<</span>EOF
        15");
        
        system("gmx editconf -f $filename\_processed.gro -o $filename\_newbox.gro -c -d 1.0 -bt cubic");
        
        system("gmx solvate -cp $filename\_newbox.gro -cs spc216.gro -o $filename\_solv.gro -p topol.top");
        
        system("gmx grompp -f ions.mdp -c $filename\_solv.gro -p topol.top -o ions.tpr");
        
        system("gmx genion -s ions.tpr -o $filename\_solv_ions.gro -p topol.top -pname NA -nname CL -neutral <span><<</span>EOF
        13");
  
        system("gmx grompp -f minim.mdp -c $filename\_solv_ions.gro -p topol.top -o em.tpr");
        
        system("gmx mdrun -v -deffnm em");
        
        system("gmx energy -f em.edr -o potential.xvg <span><<</span>EOF
        12 0");
        
        system("gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr");
        
        system("gmx mdrun -deffnm nvt");
        
        system("gmx energy -f nvt.edr -o temperature.xvg <span><<</span>EOF
        15 0");
        
        system("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr");
        
        system("gmx mdrun -deffnm npt");
        
        system("gmx energy -f npt.edr -o pressure.xvg <span><<</span>EOF
        17 0");
        
        system("gmx energy -f npt.edr -o density.xvg <span><<</span>EOF
        23");
        
        system("gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr");
        
        system("gmx mdrun -deffnm md_0_1");
        
        system("gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center <span><<</span>EOF
        1 0");
        
        system("gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns <span><<</span>EOF
        4 4");
        
        system("gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns");
    </p>
    <!-- <hr> -->
    <div class="pageName">
        <h1>ions.mdp</h1>
        <div class="divider"></div>
    </div>
    <span>
        ; ions.mdp - used as input into grompp to generate ions.tpr
        ; Parameters describing what to do, when to stop and what to save
        integrator  = steep         ; Algorithm (steep = steepest descent minimization)
        emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
        emstep      = 0.01          ; Minimization step size
        nsteps      = 50000         ; Maximum number of (minimization) steps to perform

        ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
        nstlist         = 1         ; Frequency to update the neighbor list and long range forces
        cutoff-scheme	= Verlet    ; Buffered neighbor searching 
        ns_type         = grid      ; Method to determine neighbor list (simple, grid)
        coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
        rcoulomb        = 1.0       ; Short-range electrostatic cut-off
        rvdw            = 1.0       ; Short-range Van der Waals cut-off
        pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
    </span>
    <div class="pageName">
        <h1>nvt.mdp</h1>
        <div class="divider"></div>
    </div>
    <span>
        title                   = OPLS Lysozyme NVT equilibration 
        define                  = -DPOSRES  ; position restrain the protein
        ; Run parameters
        integrator              = md        ; leap-frog integrator
        nsteps                  = 50000     ; 2 * 50000 = 100 ps
        dt                      = 0.002     ; 2 fs
        ; Output control
        nstxout                 = 500       ; save coordinates every 1.0 ps
        nstvout                 = 500       ; save velocities every 1.0 ps
        nstenergy               = 500       ; save energies every 1.0 ps
        nstlog                  = 500       ; update log file every 1.0 ps
        ; Bond parameters
        continuation            = no        ; first dynamics run
        constraint_algorithm    = lincs     ; holonomic constraints 
        constraints             = h-bonds   ; bonds involving H are constrained
        lincs_iter              = 1         ; accuracy of LINCS
        lincs_order             = 4         ; also related to accuracy
        ; Nonbonded settings 
        cutoff-scheme           = Verlet    ; Buffered neighbor searching
        ns_type                 = grid      ; search neighboring grid cells
        nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
        rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
        rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
        DispCorr                = EnerPres  ; account for cut-off vdW scheme
        ; Electrostatics
        coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order               = 4         ; cubic interpolation
        fourierspacing          = 0.16      ; grid spacing for FFT
        ; Temperature coupling is on
        tcoupl                  = V-rescale             ; modified Berendsen thermostat
        tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
        tau_t                   = 0.1     0.1           ; time constant, in ps
        ref_t                   = 300     300           ; reference temperature, one for each group, in K
        ; Pressure coupling is off
        pcoupl                  = no        ; no pressure coupling in NVT
        ; Periodic boundary conditions
        pbc                     = xyz       ; 3-D PBC
        ; Velocity generation
        gen_vel                 = yes       ; assign velocities from Maxwell distribution
        gen_temp                = 300       ; temperature for Maxwell distribution
        gen_seed                = -1        ; generate a random seed
    </span>
    <div class="pageName">
        <h1>npt.mdp</h1>
        <div class="divider"></div>
    </div>
    <span>
        title                   = OPLS Lysozyme NPT equilibration 
        define                  = -DPOSRES  ; position restrain the protein
        ; Run parameters
        integrator              = md        ; leap-frog integrator
        nsteps                  = 50000     ; 2 * 50000 = 100 ps
        dt                      = 0.002     ; 2 fs
        ; Output control
        nstxout                 = 500       ; save coordinates every 1.0 ps
        nstvout                 = 500       ; save velocities every 1.0 ps
        nstenergy               = 500       ; save energies every 1.0 ps
        nstlog                  = 500       ; update log file every 1.0 ps
        ; Bond parameters
        continuation            = yes       ; Restarting after NVT 
        constraint_algorithm    = lincs     ; holonomic constraints 
        constraints             = h-bonds   ; bonds involving H are constrained
        lincs_iter              = 1         ; accuracy of LINCS
        lincs_order             = 4         ; also related to accuracy
        ; Nonbonded settings 
        cutoff-scheme           = Verlet    ; Buffered neighbor searching
        ns_type                 = grid      ; search neighboring grid cells
        nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
        rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
        rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
        DispCorr                = EnerPres  ; account for cut-off vdW scheme
        ; Electrostatics
        coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order               = 4         ; cubic interpolation
        fourierspacing          = 0.16      ; grid spacing for FFT
        ; Temperature coupling is on
        tcoupl                  = V-rescale             ; modified Berendsen thermostat
        tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
        tau_t                   = 0.1     0.1           ; time constant, in ps
        ref_t                   = 300     300           ; reference temperature, one for each group, in K
        ; Pressure coupling is on
        pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
        pcoupltype              = isotropic             ; uniform scaling of box vectors
        tau_p                   = 2.0                   ; time constant, in ps
        ref_p                   = 1.0                   ; reference pressure, in bar
        compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
        refcoord_scaling        = com
        ; Periodic boundary conditions
        pbc                     = xyz       ; 3-D PBC
        ; Velocity generation
        gen_vel                 = no        ; Velocity generation is off 
    </span>
    <div class="pageName">
        <h1>md.mdp</h1>
        <div class="divider"></div>
    </div>
    <span>
        title                   = OPLS Lysozyme NPT equilibration 
        ; Run parameters
        integrator              = md        ; leap-frog integrator
        nsteps                  = 10000000  ; 2 * 10000000 = 20000 ps (20 ns)
        dt                      = 0.002     ; 2 fs
        ; Output control
        nstxout                 = 0         ; suppress bulky .trr file by specifying 
        nstvout                 = 0         ; 0 for output frequency of nstxout,
        nstfout                 = 0         ; nstvout, and nstfout
        nstenergy               = 5000      ; save energies every 10.0 ps
        nstlog                  = 5000      ; update log file every 10.0 ps
        nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
        compressed-x-grps       = System    ; save the whole system
        ; Bond parameters
        continuation            = yes       ; Restarting after NPT 
        constraint_algorithm    = lincs     ; holonomic constraints 
        constraints             = h-bonds   ; bonds involving H are constrained
        lincs_iter              = 1         ; accuracy of LINCS
        lincs_order             = 4         ; also related to accuracy
        ; Neighborsearching
        cutoff-scheme           = Verlet    ; Buffered neighbor searching
        ns_type                 = grid      ; search neighboring grid cells
        nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
        rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
        rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
        ; Electrostatics
        coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
        pme_order               = 4         ; cubic interpolation
        fourierspacing          = 0.16      ; grid spacing for FFT
        ; Temperature coupling is on
        tcoupl                  = V-rescale             ; modified Berendsen thermostat
        tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
        tau_t                   = 0.1     0.1           ; time constant, in ps
        ref_t                   = 300     300           ; reference temperature, one for each group, in K
        ; Pressure coupling is on
        pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
        pcoupltype              = isotropic             ; uniform scaling of box vectors
        tau_p                   = 2.0                   ; time constant, in ps
        ref_p                   = 1.0                   ; reference pressure, in bar
        compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
        ; Periodic boundary conditions
        pbc                     = xyz       ; 3-D PBC
        ; Dispersion correction
        DispCorr                = EnerPres  ; account for cut-off vdW scheme
        ; Velocity generation
        gen_vel                 = no        ; Velocity generation is off </span>
</body>
</html>
<style>
    .pageName h1{
      font-size:2.5rem;
   }
   .pageName {
    width: 100%;
    border-bottom-style:solid; 
    border-bottom-width:medium;
    margin-bottom: 5vh;
    margin-top:0;
    padding-left: 5%;
}
.pageName {
        padding-left: 0%;
        margin-top:2vh;
    }
    span{
        white-space:pre-wrap;
    }
</style>
