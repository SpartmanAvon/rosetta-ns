import random
import rosetta
from toolbox import pose_from_rcsb
from heapq import nsmallest
import cPickle
import math
import time
from collections import deque


class MCNS:

    def __init__(self, pdb, centroid=False, pdb_file='', frag=False, nine_mer=False, local=False, local_size=3,
                 full=False, rosetta_refinement=False):
        """ :param pdb: :type string: pdb ID of the protein to be folded
            :param centroid: :type boolean: Option for use of centroid model
        """
        self.loops = 0                                    # Stores generation for which energy score was last calculated
        self.scores = {}                                  # Dictionary container for current gen genomes/scores
        self.scores_list = []                             # List container of current gen scores for search
        self.gen_added = 0                                # Last gen in which a point was added to novelty archive
        self.threshold = 10                               # Novelty threshold for which point is added to archive
        self.acceptance_threshold = 100                   # Novelty threshold for which move is accepted automatically
        self.num_added = 0                                # Number of points added to novelty archive
        self.switch = False                               # All atom switch
        self.temperature = 5                              # Monte Carlo temperature
        self.mover_range = 10                             # +-range of the angle in degrees in which mover moves residue
        self.local_size = local_size                      # For local mover, size of fragment to move
        self.local = local                                # Whether to use local mover
        self.novelty_archive = deque()                    # Initialize novelty archive
        self.centroid = centroid                          # If true use centroid scoring
        self.last_lowest = 0                              # For use in novelty loop
        self.last_lowest_10 = 0                           # For use in clear main loop
        self.frag = frag                                  # If true use frag mover
        self.rosetta_refinement = rosetta_refinement      # If true refine rosetta fold

        # Rosetta inits
        rosetta.init()                                    # Initialize rosetta libraries
        pose_native = pose_from_rcsb(pdb)                 # Create rosetta pose of natively folded protein from pdb file
        sequence = pose_native.sequence()                 # Get sequence of protein
        self.scorefxn = rosetta.get_fa_scorefxn()         # Create the rosetta energy score function for all atom
        
        if pdb_file != '':
            self.pose = rosetta.pose_from_pdb(pdb_file)   # If a starting pdb is given search from this pose
        elif rosetta_refinement:                          # If rosetta refinement, start from fastrelax structure
            self.pose = rosetta.pose_from_sequence(sequence)
            relax = rosetta.FastRelax()
            relax.set_scorefxn(self.scorefxn)
            relax.apply(self.pose)
        else:
            self.pose = rosetta.pose_from_sequence(sequence)  # Create the rosetta pose that will be manipulated
            
        if centroid:                                      # Switch pose to centroid if centroid option is true
            switch = rosetta.SwitchResidueTypeSetMover("centroid")
            switch.apply(self.pose)
        self.c_size = len(sequence)*2                     # Number of residues * 2 (phi and psi for each residue)
        self.native_energy = self.scorefxn(pose_native)   # Energy of the natively folded protein
        
        if centroid:                                      # Switch rosetta score function if centroid
            self.scorefxn = rosetta.create_score_function('score3')
        self.conformation = []
        
        i = 1
        while i <= len(sequence):
            self.conformation.append(self.pose.phi(i))
            self.conformation.append(self.pose.psi(i))
            i += 1

        self.mc_energy = self.scorefxn(self.pose) + 500   # Energy to be used as minimal criteria
        self.lowest = self.scorefxn(self.pose)            # Lowest energy in archive
        
        if frag:
            if nine_mer:
                fragset = rosetta.ConstantLengthFragSet(9)
                fragset.read_fragment_file("aat000_09_05-1.200_v1_3")
            else:
                fragset = rosetta.ConstantLengthFragSet(3)
                fragset.read_fragment_file("aat000_03_05-1.200_v1_3")
            movemap = rosetta.MoveMap()
            movemap.set_bb(True)
            self.mover_3mer = rosetta.ClassicFragmentMover(fragset, movemap)

        if local:                                         # For local, initialize na with appropriate number of deques
            self.novelty_archive = [deque() for i in range(self.c_size/2/self.local_size)]

        self.full = full                                  # If true use full mover

    def frag_mover(self):

        pose_last = self.pose

        # Apply fragment mover
        self.mover_3mer.apply(self.pose)

        score = self.scorefxn(self.pose)

        novelty = 0

        # print str(score)

        if score > self.mc_energy:
            novelty = 0
        else:
            k_nearest = nsmallest(15, self.novelty_archive, key=lambda l: abs(l-score))
            for kn in k_nearest:
                novelty += abs(score - kn)
            if len(self.novelty_archive) > 15:
                novelty /= 15
            elif len(self.novelty_archive) != 0:
                novelty /= len(self.novelty_archive)
            elif len(self.novelty_archive) == 0:
                self.gen_added = self.loops
                self.novelty_archive.append(score)
                self.num_added += 1

        # print str(novelty)
        # print str(clash)

        # Decide whether to accept move and add to novelty archive
        if novelty >= self.threshold:
            self.novelty_archive.append(score)
            self.novelty_archive = sorted(self.novelty_archive)
            
            # Save to pdb if lowest
            if score == self.novelty_archive[0]:
                self.pose.dump_pdb("lowest.pdb")
                self.lowest = score
                self.pose.dump_pdb("lowest_backup.pdb")
            self.gen_added = self.loops
            self.num_added += 1
            if len(self.novelty_archive) > 10000:
                self.novelty_archive.pop()
        elif (novelty < self.acceptance_threshold) and (novelty != 0):

            # Give chance to accept move anyway
            p_acc = math.exp(-(self.threshold - novelty)/self.temperature)
            rand = random.random()
            if p_acc < rand:
                self.pose = pose_last
        elif novelty == 0:
            self.pose = pose_last

    def mover(self, n):
        c = 0
        p_s = random.random()
        p = random.randint(self.mover_range, self.mover_range)
        phi = True
        if p_s > 0.5:
            c = self.pose.phi(n)
            self.pose.set_phi(n, c + p)
        else:
            c = self.pose.psi(n)
            self.pose.set_psi(n, c + p)
            phi = False

        # print "Calculating energy..."

        score = self.scorefxn(self.pose)

        # print "Done \n"

        # print "Evaluating novelty... "
        novelty = 0
        clash = False

        # Calculate novelty based on k-nearest neighbors
        if score > self.mc_energy:
            novelty = 0
        else:

            # Check for large energy indicative of a clash. If found, set novelty score to 0
            if self.pose.energies().residue_total_energies(n)[rosetta.core.scoring.total_score] > 50000:
                clash = True
            if not clash:
                k_nearest = nsmallest(15, self.novelty_archive, key=lambda l: abs(l-score))
                for kn in k_nearest:
                    novelty += abs(score - kn)
                if len(self.novelty_archive) > 15:
                    novelty /= 15
                elif len(self.novelty_archive) != 0:
                    novelty /= len(self.novelty_archive)
                elif len(self.novelty_archive) == 0:
                    self.gen_added = self.loops
                    self.novelty_archive.append(score)
                    self.num_added += 1

        # print str(novelty)
        # print str(clash)

        # Decide whether to accept move and add to novelty archive
        if novelty >= self.threshold:
            self.novelty_archive.append(score)
            self.novelty_archive = sorted(self.novelty_archive)
        
            # Save to pdb if lowest
            if score == self.novelty_archive[0]:
                self.pose.dump_pdb("lowest.pdb")
                self.lowest = score
                self.pose.dump_pdb("lowest_backup.pdb")
            self.gen_added = self.loops
            self.num_added += 1
            if len(self.novelty_archive) > 10000:
                self.novelty_archive.pop()
        elif (novelty < self.acceptance_threshold) and (novelty != 0):

            # Give chance to accept move anyway
            p_acc = math.exp(-(self.threshold - novelty)/self.temperature)
            rand = random.random()
            if p_acc < rand:

                # Revert pose
                if phi:
                    self.pose.set_phi(n, c)
                else:
                    self.pose.set_psi(n, c)
        elif novelty == 0:
            if phi:
                self.pose.set_phi(n, c)
            else:
                self.pose.set_psi(n, c)
        # print "Done. \n"

    def full_mover(self):
        prev_pose = []
        n = 1
        while n <= self.c_size/2:
            c = 0
            p_s = random.random()
            p = random.randint(self.mover_range, self.mover_range)
            phi = True
            if p_s > 0.5:
                c = self.pose.phi(n)
                self.pose.set_phi(n, c + p)
            else:
                c = self.pose.psi(n)
                self.pose.set_psi(n, c + p)
                phi = False
            prev_pose.append([c, phi])

            n += 1

        # print "Calculating energy..."

        score = self.scorefxn(self.pose)

        # print "Done \n"

        # print "Evaluating novelty... "
        novelty = 0
        clash = False

        # print str(score)

        # Calculate novelty based on k-nearest neighbors
        if score > self.mc_energy:
            novelty = 0
        else:

            # Check for large energy indicative of a clash. If found, set novelty score to 0
            n = 1
            while n <= self.c_size/2:
                if self.pose.energies().residue_total_energies(n)[rosetta.core.scoring.total_score] > 50000:
                    clash = True
                n += 1
            if not clash:
                k_nearest = nsmallest(15, self.novelty_archive, key=lambda l: abs(l-score))
                for kn in k_nearest:
                    novelty += abs(score - kn)
                if len(self.novelty_archive) > 15:
                    novelty /= 15
                elif len(self.novelty_archive) != 0:
                    novelty /= len(self.novelty_archive)
                elif len(self.novelty_archive) == 0:
                    self.gen_added = self.loops
                    self.novelty_archive.append(score)
                    self.num_added += 1

        # print str(novelty)
        # print str(clash)

        # Decide whether to accept move and add to novelty archive
        if novelty >= self.threshold:
            self.novelty_archive.append(score)
            self.novelty_archive = sorted(self.novelty_archive)
        
            # Save to pdb if lowest
            if score == self.novelty_archive[0]:
                self.pose.dump_pdb("lowest.pdb")
                self.lowest = score
                self.pose.dump_pdb("lowest_backup.pdb")
            self.gen_added = self.loops
            self.num_added += 1
            if len(self.novelty_archive) > 10000:
                self.novelty_archive.pop()
        elif (novelty < self.acceptance_threshold) and (novelty != 0):

            # Give chance to accept move anyway
            p_acc = math.exp(-(self.threshold - novelty)/self.temperature)
            rand = random.random()
            if p_acc < rand:

                # Revert pose
                n = 1
                while n <= self.c_size/2:
                    res = prev_pose[n - 1]
                    if res[1]:
                        self.pose.set_phi(n, res[0])
                    else:
                        self.pose.set_psi(n, res[0])
                    n += 1
        elif novelty == 0:
            n = 1
            while n <= self.c_size/2:
                res = prev_pose[n - 1]
                if res[1]:
                    self.pose.set_phi(n, res[0])
                else:
                    self.pose.set_psi(n, res[0])
                n += 1
        # print "Done. \n"

    def local_mover(self, n, n_res):

        pose_last = self.pose
        score = 0
        clash = False

        for i in range(n_res):
            r1 = random.randint(-self.mover_range, self.mover_range)
            r2 = random.randint(-self.mover_range, self.mover_range)
            ph = self.pose.phi(n + i)
            self.pose.set_phi(n + i, ph + r1)
            ps = self.pose.psi(n + i)
            self.pose.set_psi(n + i, ps + r2)
            self.scorefxn(self.pose)
            energy = self.pose.energies().residue_total_energies(n + i)[rosetta.core.scoring.total_score]
            score += energy
            if energy > 50000:
                clash = True

        total_score = self.scorefxn(self.pose)
        if total_score < self.lowest:
            self.lowest = total_score

        novelty = 0
        novelty_index = n/self.local_size

        print str(total_score)

        if total_score > self.mc_energy:
            pass
        elif not clash:
            archive = self.novelty_archive[novelty_index]
            k_nearest = nsmallest(15, archive, key=lambda l: abs(l-score))
            for kn in k_nearest:
                novelty += abs(score - kn)
            if len(archive) > 15:
                novelty /= 15
            elif len(archive) != 0:
                novelty /= len(archive)
            elif len(archive) == 0:
                self.gen_added = self.loops
                self.novelty_archive[novelty_index].append(score)
                self.num_added += 1

        print str(novelty)
        # print str(clash)

        # Decide whether to accept move and add to novelty archive
        if novelty >= self.threshold:
            self.novelty_archive[novelty_index].append(score)
            self.novelty_archive[novelty_index] = sorted(self.novelty_archive[novelty_index])
            # Save to pdb if lowest
            if score == self.novelty_archive[novelty_index][0]:
                self.pose.dump_pdb("lowest.pdb")
                self.pose.dump_pdb("lowest_backup.pdb")
            self.gen_added = self.loops
            self.num_added += 1
            if len(self.novelty_archive[novelty_index]) > 10000:
                self.novelty_archive[novelty_index].pop()
        elif (novelty < self.acceptance_threshold) and (novelty != 0):

            # Give chance to accept move anyway
            p_acc = math.exp(-(self.threshold - novelty)/self.temperature)
            rand = random.random()
            if p_acc < rand:
                self.pose = pose_last
        elif novelty == 0:
            self.pose = pose_last



    # Novelty criterion
    def go(self):
        """ :param candidates: Current population being evaluated

            Scores the conformation by novelty.
        """

        time.clock()
        while (self.loops < 1000000) or (self.lowest > self.native_energy):

            # Do MCNS on each residue in protein
            if self.frag:
                n = 1
                while n <= self.c_size/2:
                    self.frag_mover()
                    n += 1
            elif self.local:
                n = 1
                while n < self.c_size/2/self.local_size:
                    self.local_mover(n, self.local_size)
                    n += 5
                r = self.c_size/2 % self.local_size
                self.local_mover(n, r)
            elif self.full:
                self.full_mover()
            else:
                n = 1
                while n <= self.c_size/2:
                    self.mover(n)
                    n += 1

            # Options for dynamic parameter adjustment

            # If no new points have been added in 4 loops decrease novelty threshold by 5%
            # if self.loops % 4 == 0:
            #     if self.num_added == 0:
            #         self.acceptance_threshold *= .95

            # If more than 20 points has been added in 4 loops increase novelty threshold by 20%
            # if self.loops % 4 == 0:
            #     if self.num_added > 20:
            #         self.acceptance_threshold *= 1.2

            # If no new points have been added in 1 loop increase mover range by 10
            # if self.num_added == 0:
            #     if self.mover_range < 350:
            #         self.mover_range += 10

            # If more than 1 point has been added in 1 loop decrease mover range by 5
            # if self.num_added > 1:
            #     if self.mover_range > 5:
            #         self.mover_range -= 5

            # Do a individual residue mover round if no improvement in 20 moves
            # if (self.loops % 20 == 0) and (self.loops != 0):
            #     if self.lowest == self.last_lowest_10:
            #         n = 1
            #         while n <= self.c_size/2:
            #             self.mover(n)
            #             n += 1
            #     self.last_lowest_10 = self.lowest

            # Print lowest e_score in archive compared to native
            if len(self.novelty_archive) > 0:
                print "Lowest energy: " + str(self.lowest) + "\nNative energy: " + str(self.native_energy) + \
                    "\nNovelty points added: " + str(self.num_added) + "\n"

                # If lowest e_score is lower then 10000 convert to all-atom, switch energy function, and change minimal
                # criteria to new all-atom energy (for centroid)
                if (self.lowest < 10000) and not self.switch and self.centroid:
                    # self.threshold = 10
                    switch_type = rosetta.SwitchResidueTypeSetMover("fa_standard")
                    switch_type.apply(self.pose)
                    self.scorefxn = rosetta.get_fa_scorefxn()
                    self.mc_energy = self.scorefxn(self.pose)
                    self.switch = True

            # Clear novelty archive after 10 loops and set mc energy to lowest if no improvement
            # if (self.loops % 10 == 0) and (self.loops != 0):
            #     if self.lowest == self.last_lowest_10:
            #         print("Reloading...")
            #         self.pose = rosetta.pose_from_pdb('lowest.pdb')
            #         self.novelty_archive = deque()
            #         self.mc_energy = self.scorefxn(self.pose) + 500
            #         self.acceptance_threshold = 100
            #     self.last_lowest_10 = self.lowest

            # Decrease temp by 5 if no progress after 20 loops
            # if self.loops % 10 == 0:
            #     if self.lowest == self.last_lowest:
            #         self.temperature += 5
            #     self.last_lowest = self.lowest

            self.loops += 1
            self.num_added = 0
            print str(self.loops) + " iterations."
            # print "Threshold: " + str(self.threshold)
        print ("Time elapsed: " + str(time.clock()))

    def run(self):
        """ Run evolution on the protein
        """
        # Dump novelty archive to pickle
        output = open("confs.pkl", "wb")
        cPickle.dump(self.novelty_archive, output)
        output.close()

        # Print lowest energy obtained
        ns = sorted(self.novelty_archive)
        print(str(ns[0][0]))
