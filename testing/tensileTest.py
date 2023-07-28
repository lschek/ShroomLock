import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "cmss",
    "font.size": 20,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 20
})

class experiment:
    """
    This class is for every combination (e.g. same parts). Displacement and force will
    be stored in an array that contains results for all takes of that experiment.
    """
    def __init__(self, title):
        self.title = title

    def read_data(self, n_takes=1):
        """
        Reads data from .xlsx files in data folder.
        """
        self.n_takes = n_takes
        # TODO add assert to check if file exists
        self.displacement = np.zeros(n_takes, dtype=object)
        self.force = np.zeros(n_takes, dtype=object)
        for n in range(n_takes):
            number = str(n + 1).zfill(2) # leading zeros for all number < 10
            df = pd.read_excel("../data/" + self.title + "_" + number + ".xlsx", sheet_name=2, skiprows=1)

            self.displacement[n] = np.array(df["Standardweg"][1:], dtype=float)
            self.force[n] = np.array(df["Standardkraft"][1:], dtype=float)

    def interpolate_data(self, n_interpPoints, zmax):
        self.interpolated_displacement = np.linspace(0, zmax, n_interpPoints)
        mean_forces = np.zeros(self.n_takes, dtype=object)
        for n in range(1, self.n_takes):
            mean_forces[n] = np.interp(self.interpolated_displacement, self.displacement[n], self.force[n])

        self.mean_force = np.mean(mean_forces[1:], axis=0)
        self.std_force = np.std(mean_forces[1:], axis=0)
    
    def make_plot(self, xmax, ymax):
        self.fig, axs = plt.subplots(1,2,figsize=(16,8))

        for n in range(self.n_takes):
            axs[1].plot(self.displacement[n], self.force[n], label=f"Measurement {n + 1}")
            #axs[1].plot(self.displacement[n], self.force[n], label=f"Messung {n + 1}")

        
        meancolor = "lightskyblue"
        axs[0].plot(self.interpolated_displacement, self.mean_force, color=meancolor, label=f"Mean of all measurements")
        #axs[0].plot(self.interpolated_displacement, self.mean_force, color=meancolor, label=f"Mittelung aller Messungen")
        axs[0].fill_between(self.interpolated_displacement, self.mean_force-self.std_force, self.mean_force + self.std_force, color=meancolor, alpha=0.3)

        # Turn on the grid
        for ax in axs:
            ax.grid(True, linestyle='dotted')
            ax.set_xlabel(r"Displacement $u_z$ in mm")
            ax.set_ylabel(r"Force $F$ in N")
            ax.set_xlim(0, xmax) #1.02 * self.interpolated_displacement[-1])
            ax.set_ylim(0, ymax) #1.02 * np.max([np.max(arr) for arr in self.force]))
            ax.legend(loc='upper right')

        self.fig.suptitle(self.title)
        plt.tight_layout()

class mechanism:
    def __init__(self, n_experiments):
        self.n_experiments = n_experiments
        self.experiments = [experiment(f"TPU_V7_M{n + 1}") for n in range(n_experiments)]
        
        for experiment_n in self.experiments:
            experiment_n.read_data(4) # TODO generalize!
        self.max_displacement = np.max([np.max([np.max(arr) for arr in experiment_n.displacement]) for experiment_n in self.experiments])
        self.max_force = np.max([np.max([np.max(arr) for arr in experiment_n.force]) for experiment_n in self.experiments])

        for experiment_n in self.experiments:
            experiment_n.interpolate_data(n_interpPoints=30, zmax=self.max_displacement)
            experiment_n.make_plot(xmax=1.1*self.max_displacement, ymax=1.1*self.max_force)
            experiment_n.fig.savefig(f"{experiment_n.title}-results.pdf", format="pdf")
            experiment_n.fig.savefig(f"{experiment_n.title}-results.svg")

    def make_plot(self, xmax=None, ymax=None):
        if xmax == None:
            xmax = 1.1 * self.max_displacement
        if ymax == None:
            ymax = 1.1 * self.max_force

        self.mean_force = np.mean([experiment_n.mean_force for experiment_n in self.experiments], axis=0)
        self.std_force = np.std([experiment_n.mean_force for experiment_n in self.experiments], axis=0)
        self.displacement = self.experiments[0].interpolated_displacement # TODO

        self.fig, axs = plt.subplots(1,2,figsize=(16,8))

        meancolor = "lightskyblue" #TODO
        axs[0].plot(self.displacement, self.mean_force, color=meancolor, label="Mean of means of all mechanisms")
        #axs[0].plot(self.displacement, self.mean_force, color=meancolor, label="Mittelung aller Mechanismen")
        axs[0].fill_between(self.displacement, self.mean_force-self.std_force, self.mean_force + self.std_force, color=meancolor, alpha=0.3)
        #u_foo = np.array([0.16, 0.41, 0.74, 1.12, 1.48, 1.81, 2.11, 2.38, 2.65, 2.94])

        #F_foo = 4 * np.array([0.5, 1.,  1.5, 2.,  2.5, 3.,  3.5, 4.,  4.5, 5. ])
        #axs[0].plot(u_foo, F_foo, label="Numerical analysis, without slip")
        for n, experiment_n in enumerate(self.experiments):
            axs[1].plot(experiment_n.interpolated_displacement, experiment_n.mean_force, label=f"Mean of mechanism {n + 1}")
            #axs[1].plot(experiment_n.interpolated_displacement, experiment_n.mean_force, label=f"Mechanismus {n + 1}")
            axs[1].fill_between(experiment_n.interpolated_displacement, experiment_n.mean_force-experiment_n.std_force, experiment_n.mean_force + experiment_n.std_force, alpha=0.3)

                # Turn on the grid
        for ax in axs:
            ax.grid(True, linestyle='dotted')
            ax.set_xlabel(r"Displacement $u_z$ in mm")
            ax.set_ylabel(r"Force $F$ in N")
            ax.set_xlim(0, xmax)
            ax.set_ylim(0, ymax)
            #ax.legend(loc='upper right')
        axs[0].legend(loc='upper right')
        self.fig.suptitle("Mean result from all mechanisms")
        plt.tight_layout()


        
if __name__=="__main__":
    # experiment1 = experiment("TPU_V7_M1")
    # experiment1.read_data(4)
    # experiment1.make_plot()
    # experiment1.fig.savefig("output.pdf", format="pdf")

    # for foo in range(4):
    #     title = f"TPU_V7_M{foo + 1}"
    #     experiment_foo = experiment(title)
    #     experiment_foo.read_data(4)
    #     experiment_foo.interpolate_data(30) # TODO all data should be interpolated over the same x interval as we are later calculating the mean
    #     experiment_foo.make_plot(5.5, 18.5) # TODO When it comes to plotting, keep in mind the biggest numbers!
    #     experiment_foo.fig.savefig(f"{title}-results.pdf", format="pdf")

    foo = mechanism(8)
    foo.make_plot()
    foo.fig.savefig("TPU_V7-results.pdf", format="pdf")
    foo.fig.savefig("TPU_V7-results.svg")