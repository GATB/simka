
import os, sys, shutil, glob
os.chdir(os.path.split(os.path.realpath(__file__))[0])

input_filename = sys.argv[1]
nb_boostraps = int(sys.argv[2])
output_dir_temp = os.path.join(sys.argv[3], "__temp__")
SEQUENCING_EFFORT_MAX = 1000000
SEQUENCING_EFFORT = [10000, 100000]#, 10000000, 100000000]






simka_tmp_dir = os.path.join(output_dir_temp, "simka_temp")
boostrap_results_dir = os.path.join(output_dir_temp, "boostrap_results")
r_input_dir = os.path.join(output_dir_temp, "r_input_dir")
r_result_dir = os.path.join(output_dir_temp, "result_figures")




simka_out_dir = os.path.join(output_dir_temp, "simka_results")
simka_command = "../../build/bin/simka "
simka_command += " -in " + input_filename
simka_command += " -out-tmp " + simka_tmp_dir
simka_command += " -kmer-size 21 "
simka_command += " -abundance-min 0 "
simka_command += " -subsampling-kind 1 "



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
        filename = os.path.join(output_dir_temp, "simka_subsampling_setup.txt")
        command = "../../build/bin/simka -in " + input_filename + " -out-tmp " + output_dir_temp + " -subsampling-setup"
        command += " > " + filename
        #print command
        #exit(1)
        #os.system(command)

        #for line in open(filename, "r"):
        #    if "Reference dataset ID" in line:
        #        self.subsampling_reference_dataset_ID = int(line.strip().replace(" ", "").replace("ReferencedatasetID:", ""))
        #   if "Subsampling space" in line:
        #       self.subsampling_max_reads = int(line.strip().replace(" ", "").replace("Subsamplingspace(reads):", ""))
        #       #break

        #print("Reference dataset ID: " + str(self.subsampling_reference_dataset_ID))
        #print("Subsampling space: " + str(self.subsampling_max_reads))

    def compute_truth(self):

        command = simka_command
        command += " -subsampling-space " + str(SEQUENCING_EFFORT_MAX)
        command += " -subsampling-ref-id 0 "
        command += " -subsampling-nb-reads 0 "
        #command += " -max-reads " + str(self.subsampling_max_reads)

        output_dir = os.path.join(output_dir_temp, "truth_results")

        self.run_simka(command, output_dir)


    def subsample(self):

        for sequencing_effort in SEQUENCING_EFFORT:

            #nb_reads_to_pick = int((self.subsampling_kmer_space * percent) / float(100))

            command = simka_command
            command += " -subsampling-space " + str(SEQUENCING_EFFORT_MAX)
            command += " -subsampling-ref-id 0 "
            command += " -subsampling-nb-reads " + str(sequencing_effort)

            for i in range(0, nb_boostraps):
                boostrap_out_dir = os.path.join(boostrap_results_dir, "pass_" + str(sequencing_effort) + "_" + str(i))

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

s = ComputeBootstraps()
s.execute()
