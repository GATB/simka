import os, sys, shutil, glob
os.chdir(os.path.split(os.path.realpath(__file__))[0])

input_filename = sys.argv[1]
nb_boostraps = int(sys.argv[2])
subsampling_space = int(sys.argv[4])
output_dir_temp = os.path.join(sys.argv[3], str(subsampling_space), "__temp__")

IS_PICK_PERCENTAGE = True
NB_READS_TO_PICK = [0.1, 0.5, 1, 10, 20, 30, 40, 50, 80]






simka_tmp_dir = os.path.join(output_dir_temp, "simka_temp")
boostrap_results_dir = os.path.join(output_dir_temp, "boostrap_results")
r_input_dir = os.path.join(output_dir_temp, "r_input_dir")
r_result_dir = os.path.join(output_dir_temp, "result_figures")




simka_out_dir = os.path.join(output_dir_temp, "simka_results")
simka_command = "../../build/bin/simka "
simka_command += " -in " + input_filename
simka_command += " -out-tmp " + simka_tmp_dir
simka_command += " -kmer-size 31 "
simka_command += " -abundance-min 0 "
simka_command += " -max-memory 10000 "
simka_command += " -nb-cores 8 "
simka_command += " -subsampling-kind 1 "
#simka_command += " -simple-dist "
#simka_command += " -complex-dist "

def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)



create_dir(output_dir_temp)
create_dir(boostrap_results_dir)
create_dir(r_input_dir)
create_dir(r_result_dir)


class ComputeBootstraps():

    def execute(self):
        self.setup()
        self.compute_truth()
        self.subsample()

    def setup(self):
        #filename = os.path.join(output_dir_temp, "simka_subsampling_setup.txt")
        #command = "../../build/bin/simka -in " + input_filename + " -out-tmp " + output_dir_temp + " -subsampling-setup"
        #command += " > " + filename
        #print command
        #exit(1)
        #os.system(command)

        #for line in open(filename, "r"):
        #    if "Reference dataset ID" in line:
        #        self.subsampling_reference_dataset_ID = int(line.strip().replace(" ", "").replace("ReferencedatasetID:", ""))
        #    if "Subsampling space" in line:
        #        self.subsampling_max_reads = int(line.strip().replace(" ", "").replace("Subsamplingspace(reads):", ""))
        #        #break

        self.subsampling_reference_dataset_ID = 0
        self.subsampling_max_reads = subsampling_space

        print("Reference dataset ID: " + str(self.subsampling_reference_dataset_ID))
        print("Subsampling space: " + str(self.subsampling_max_reads))

    def compute_truth(self):

        #command = simka_command
        #command += " -subsampling-space " + str(self.subsampling_kmer_space)
        #command += " -max-reads " + str(self.subsampling_max_reads)

        command = simka_command
        command += " -subsampling-space " + str(self.subsampling_max_reads)
        command += " -subsampling-ref-id " + str(self.subsampling_reference_dataset_ID)
        command += " -subsampling-nb-reads " + str(self.subsampling_max_reads)

        output_dir = os.path.join(output_dir_temp, "truth_results")

        self.run_simka(command, output_dir)


    def subsample(self):

        for nb_reads_to_pick in NB_READS_TO_PICK:

            percent_nb_reads = nb_reads_to_pick
            nb_reads_to_pick = self.get_nb_reads_to_pick(nb_reads_to_pick)

            if nb_reads_to_pick >= self.subsampling_max_reads: continue

            #nb_reads_to_pick = int((self.subsampling_kmer_space * percent) / float(100))

            command = simka_command
            command += " -subsampling-space " + str(self.subsampling_max_reads)
            command += " -subsampling-ref-id " + str(self.subsampling_reference_dataset_ID)
            command += " -subsampling-nb-reads " + str(nb_reads_to_pick)

            for i in range(0, nb_boostraps):

                if IS_PICK_PERCENTAGE:
                    boostrap_out_dir = os.path.join(boostrap_results_dir, "pass_" + str(percent_nb_reads) + "_" + str(i))
                else:
                    boostrap_out_dir = os.path.join(boostrap_results_dir, "pass_" + str(nb_reads_to_pick) + "_" + str(i))

                self.run_simka(command, boostrap_out_dir)
                #os.system(command + " > " + os.path.join(boostrap_out_dir, "log.txt"))

                #print boostrap_out_dir
                #shutil.move(simka_out_dir, boostrap_out_dir)
                #exit(0)

    def run_simka(self, command, output_dir):
        if self.is_pass_already_computed(output_dir): return

        #print output_dir
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)

        if os.path.exists(simka_tmp_dir):
            shutil.rmtree(simka_tmp_dir)

        command += " -out " + output_dir

        print command + " > " + os.path.join(output_dir, "log.txt")
        os.system(command + " > " + os.path.join(output_dir, "log.txt"))
        #shutil.move(simka_out_dir, output_dir)


    def is_pass_already_computed(self, simka_results_dir):
        return os.path.exists(simka_results_dir) and len(glob.glob(os.path.join(simka_results_dir, "mat_*"))) > 0

    def get_nb_reads_to_pick(self, nb_reads_to_pick):
        if IS_PICK_PERCENTAGE:
            return (nb_reads_to_pick * self.subsampling_max_reads) / 100
        else:
            return nb_reads_to_pick

class ComputeBootstrapsStats():

    def __init__(self):
        pass

    def execute(self):
        data = {}

        for dir in glob.glob(os.path.join(boostrap_results_dir, "pass_*")):
            if "asym" in dir: continue

            #print dir
            basename = os.path.basename(dir)
            dummy, percent, passID = basename.split("_")
            #print percent

            matrix_filenames = glob.glob(os.path.join(dir, "mat_*"))
            for matrix_filename in matrix_filenames:
                fields = os.path.basename(matrix_filename).split(".")[0].split("_")
                distance_name = fields[1] + "_" + fields[2]

                if not distance_name in data:
                    data[distance_name] = {}

                if not percent in data[distance_name]:
                    data[distance_name][percent] = []

                data[distance_name][percent].append(matrix_filename)

        #print data
        distance_names = []
        distance_filenames = glob.glob(os.path.join(output_dir_temp, "truth_results", "mat_*"))
        for distance_filename in distance_filenames:
            if "asym" in distance_filename: continue
            distance_names.append(os.path.basename(distance_filename).split(".")[0].replace("mat_", ""))

        input_filename_R = os.path.join(r_input_dir, "input_bootstrap.txt")
        #print data["abundance_braycurtis"]
        #distance_name = "abundance_braycurtis"
        for distance_name in distance_names:
            create_dir(os.path.join(r_result_dir, distance_name))

            input_R_file = open(input_filename_R, "w")

            for percent, matrix_filenames  in data[distance_name].items():
                input_R_file.write(percent)
                for filename in matrix_filenames:
                    input_R_file.write(" " + filename)
                input_R_file.write("\n")
                #print percent, matrix_filenames
                #print matrix_filenames
                #data[percent].append()
            input_R_file.close()

            os.system("Rscript subsampling_stats.r " + distance_name + " " + input_filename_R + " " + output_dir_temp)

s = ComputeBootstraps()
s.execute()

#s = ComputeBootstrapsStats()
#s.execute()
