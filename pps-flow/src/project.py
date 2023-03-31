"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os


class MyProject(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="gpu",
            help="Specify the partition to submit to."
        )


class R2(DefaultSlurmEnvironment):
    hostname_pattern = "r2"
    template = "r2.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="shortgpuq",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )

# Definition of project-related labels (classification)
@MyProject.label
def sampled(job):
    return job.doc.get("done")


@MyProject.label
def initialized(job):
    pass


@directives(executable="python -u")
@directives(ngpu=1)
@MyProject.operation
@MyProject.post(sampled)
def sample(job):
    import hoomd_polymers
    from hoomd_polymers.systems import Pack
    import hoomd_polymers.forcefields
    from hoomd_polymers.forcefields import OPLS_AA_PPS
    from hoomd_polymers.sim import Simulation
    from hoomd_polymers.molecules import PPS
    from cmeutils.sampling import is_equilibrated, equil_sample 

    with job:
        print("JOB ID NUMBER:")
        print(job.id)

        system = Pack(
                molecule=PPS,
                density=job.sp.density,
                n_mols=job.sp.n_compounds,
                mol_kwargs = {
                    "length": job.sp.chain_lengths,
                },
                packing_expand_factor=5
        )

        system.apply_forcefield(
                forcefield=OPLS_AA_PPS(),
                make_charge_neutral=True,
                remove_charges=job.sp.remove_charges,
                remove_hydrogens=job.sp.remove_hydrogens
        )

        job.doc.ref_distance = system.reference_distance
        job.doc.ref_mass = system.reference_mass
        job.doc.ref_energy = system.reference_energy

        gsd_path = os.path.join(job.ws, "trajectory.gsd")
        log_path = os.path.join(job.ws, "sim_data.txt")

        sim = Simulation(
            initial_state=system.hoomd_snapshot,
            forcefield=system.hoomd_forcefield,
            dt=job.sp.dt,
            gsd_write_freq=job.sp.gsd_write_freq,
            gsd_file_name=gsd_path,
            log_file_name=log_path,
            log_write_freq=job.sp.log_write_freq
        )
        sim.pickle_forcefield(job.fn("pps_forcefield.pickle"))

        target_box = system.target_box*10/job.doc.ref_distance
        job.doc.target_box = target_box
        job.doc.real_timestep = sim.real_timestep.to("fs")
        job.doc.real_timestep_units = "fs"
        print("-----------------")
        print("Running shrink step...")
        print("-----------------")
        sim.run_update_volume(
                final_box_lengths=target_box,
                n_steps=job.sp.shrink_steps,
                period=job.sp.shrink_period,
                tau_kt=job.sp.tau_kt,
                kT=job.sp.shrink_kT
        )
        print("-----------------")
        print("Shrink step finished; running NVT...")
        print("-----------------")
        sim.run_NVT(kT=job.sp.kT, n_steps=job.sp.n_steps, tau_kt=job.sp.tau_kt)

        job.doc.shrink_cut = int(job.sp.shrink_steps/job.sp.log_write_freq) 
        run_num = 1
        equilibrated = False
        while not equilibrated:
        # Open up log file, see if pressure is equilibrated
            data = np.genfromtxt(job.fn("sim_data.txt"), names=True)
            pressure = data["mdcomputeThermodynamicQuantitiespressure"]
            equilibrated = is_equil(pressure[job.doc.shrink_cut:], threshold_neff=5000)
            print("-----------------")
            print(f"Not yet equilibrated. Starting run {run_num + 1}.")
            print("-----------------")
            sim.run_NVT(kT=job.sp.kT, n_steps=job.sp.extra_steps, tau_kt=job.sp.tau_kt)
            run_num += 1

        print("-----------------")
        print("Is equilibrated; starting sampling...")
        print("-----------------")
        data = np.genfromtxt(job.fn("sim_data.txt"), names=True)
        pressure = data["mdcomputeThermodynamicQuantitiespressure"]
        uncorr_sample, uncorr_indices, prod_start, ineff, Neff = equil_sample(
                pressure[job.doc.shrink_cut:],
                threshold_fraction=0.50,
                threshold_neff=5000
        )
        job.doc.uncorr_indices = list(uncorr_indices)
        job.doc.prod_start = prod_start
        job.doc.average_pressure = np.mean(uncorr_sample)
        job.doc.pressure_std = np.std(uncorr_sample)
        job.doc.pressure_sem = np.std(uncorr_sample) / (len(uncorr_sample)**0.5)

        job.doc.total_steps = job.sp.n_steps + (run_num * job.sp.extra_steps)
        job.doc.total_time_ns = job.doc.total_steps*job.doc.real_timestep.to("ns")
        job.doc.n_runs = run_num
        job.doc.done = True


if __name__ == "__main__":
    MyProject().main()
