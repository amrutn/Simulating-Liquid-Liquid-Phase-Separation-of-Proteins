import numpy as np
import tqdm
import time
from helper import closestDistanceBetweenLines

""" All units are in um, uM and s. diffusion_var_per_sec is in um^2/s. Each lattice point is the size of a molecule"""



class Aggregate:
    """
    Represents an aggregate of one or more molecules.
    """

    molecular_radius = 0.01 #um
    def __init__(self, num_molecules, coordinates, diffusive_const, is_droplet = False):
        """
        Initializes a molecular aggregate.

        Parameters
        ----------
        num_molecules : int
            Number of molecules in the aggregate
        radius : float
            The radius of the aggregate
        coordinates : array of three floats
            Represents the cartesian location of the center of mass of the aggregate
        diffusive_const : float
            A measure of how difficult it is for an aggregate to move in the surrounding medium.
        is_droplet : bool
            True if the aggregate is large enough to be considered a droplet

        """
        self.num_molecules = num_molecules
        self.radius = (num_molecules)**(1/3) * self.molecular_radius
        self.coordinates = coordinates
        self.diffusion_var_per_sec = 6 * diffusive_const/self.radius #See solution to diffusion equation and its variance
        self.diffuse_time = 0

    def volume(self):
        """
        Returns volume of the aggregate
        """
        return 4/3 * np.pi * self.radius**3

    def coords(self):
        """
        Returns coordinates of the aggregate.
        """
        return np.array(self.coordinates)
    
    def is_droplet(self):
        """
        Returns whether or not the aggregate is large enough to be considered a droplet.
        """
        return (self.volume() > 0.0001)

    def update_diffuse_time(self):
        """
        Increases the amount of time without diffusion by 1 for the aggregate.
        """
        self.diffuse_time += 1

    def update_coords(self):
        """
        Moves the aggregate in a random direction by a distance set by the diffusion constant
        """
        delta = np.random.normal(0, np.sqrt(self.diffusion_var_per_sec), 3)
        self.coordinates = np.add(self.coords(), delta)
        return self.coords()



        
class Sample:

    """
    A class that represents a sample of dissolved TRF1 molecules. 
    Simulates the movement and aggregation of these molecules 
    into droplets. 
    """

    def __init__(self, sample_volume, diffusive_const, concentration):

        """
        Initializes the sample.
        """
        self.diffusive_const = diffusive_const
        self.num_molecules = concentration * 6 * 10**2 * sample_volume
        self.sample_volume = sample_volume
        self.aggregates = []
        i = 0
        while i < self.num_molecules:
            coords = sample_volume**(1/3) * np.random.rand(3)
            self.aggregates.append(Aggregate(1, coords, diffusive_const))
            i += 1
    def num_droplets(self):
        """
        Returns the number of aggregates in self that are large enough to be considered droplets.
        """
        num_droplets = 0
        for agg in self.aggregates:
            if agg.is_droplet():
                num_droplets += 1
        return num_droplets
    
    def simulate(self, timesteps):
        """
        Simulates the placement, movement and merging of individual molecules

        Which aggregates merge is determined by simulating the movement of the aggregates
        in random directions and checking whether or not they intersect. 
    
        """
        for i in tqdm.tqdm(range(timesteps)):
            num_aggs = len(self.aggregates)
            path_list = []
            for agg in self.aggregates:
                prev_coords = agg.coords()
                agg.update_diffuse_time()
                coords = agg.update_coords()
                path_list.append([prev_coords, coords])
            i = 0
            while i < num_aggs:
                j = i + 1
                while j < num_aggs:

                    pA, pB, path_distance = closestDistanceBetweenLines(path_list[i][0], path_list[i][1], path_list[j][0], path_list[j][1])
                    if path_distance <= self.aggregates[i].radius + self.aggregates[j].radius:
                        num_merged_molecules = self.aggregates[i].num_molecules + self.aggregates[j].num_molecules
                        new_coords = np.add(pA, pB)/2
                        merged = Aggregate(num_merged_molecules, new_coords, self.diffusive_const)
                        self.aggregates.append(merged)
                        self.aggregates.pop(j)
                        self.aggregates.pop(i)
                        path_list.pop(j)
                        path_list.pop(i)
                        num_aggs -= 2
                        j = i
                    j += 1
                i += 1
        return self.aggregates


