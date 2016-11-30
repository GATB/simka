
import os, re


class SimkaDatabase():

	def __init__(self, dirname):
		self.dirname = dirname

		#Keep order of entries
		self.entries = []

		#Hold information of entries
		self.entries_infos = {}

		self.create_dirs()


		self.load_settings()
		self.load()
		#print(self.items)

	def load(self):
		#self.database_filename_entries = os.path.join(self.dirname, "simka_database_entries.bin")
		self.database_filename = os.path.join(self.dirname, "simka_database.csv")

		if os.path.exists(self.database_filename):
			self.database_file = open(self.database_filename, "r")
			self.database_file.next() #skip header
			for line in self.database_file:
				line = line.strip()
				if line == "": continue
				id, relative_filename = line.split(";")
				self.entries.append(id)
				self.entries_infos[id] = relative_filename
			self.database_file.close()
		else:

			self.database_file = open(self.database_filename, "w")
			self.database_file.write("ID;Kmer_Spectrum_Filename\n")
			self.database_file.close()


	def load_settings(self):
		settings_filename = os.path.join(self.dirname, "settings.txt")
		settings_file = open(settings_filename, "r")
		lines = []
		for line in settings_file:
			line = line.strip()
			if line == "": continue
			lines.append(line)
		self._kmerSize = int(lines[0].replace(" ", "").replace("kmer-size:", ""))
		self._nbPartitions = int(lines[1].replace(" ", "").replace("nb-partitions:", ""))
		self._abundanceMin = int(lines[2].replace(" ", "").replace("abundance-min:", ""))
		self._abundanceMax = int(lines[3].replace(" ", "").replace("abundance-max:", ""))
		self._maxReads = int(lines[4].replace(" ", "").replace("max-reads:", ""))
		self._minReadSize = int(lines[5].replace(" ", "").replace("min-read-size:", ""))
		self._minShannonIndex = int(lines[6].replace(" ", "").replace("min-shannon-index:", ""))
		self._computeSimpleDist = (True if (lines[7].replace(" ", "").replace("simple-dist:", "") == "1") else False)
		self._computeComplexDist = (True if (lines[8].replace(" ", "").replace("complex-dist:", "") == "1") else False)

	def create_dirs(self):
		self.kmer_spectrums_relative_dir = "kmer_spectrums"
		#self.distance_computation_dir = "distance_computation"

		kmer_spectrum_absolute_dir = os.path.join(self.dirname, self.kmer_spectrums_relative_dir)
		if not os.path.exists(kmer_spectrum_absolute_dir):
			os.mkdir(kmer_spectrum_absolute_dir)


		self.mergeKmerSpectrumRelativeDir = os.path.join(self.kmer_spectrums_relative_dir, "__merged__")
		self.mergeKmerSpectrumAbsDir = os.path.join(self.dirname, self.mergeKmerSpectrumRelativeDir)
		if not os.path.exists(self.mergeKmerSpectrumAbsDir): os.mkdir(self.mergeKmerSpectrumAbsDir)

		dir = os.path.join(self.dirname, "distance")
		if not os.path.exists(dir): os.mkdir(dir)

		dir = os.path.join(self.dirname, "distance", "temp_parts")
		if not os.path.exists(dir): os.mkdir(dir)

		dir = os.path.join(self.dirname, "distance", "matrix_binary")
		if not os.path.exists(dir): os.mkdir(dir)

		dir = os.path.join(self.dirname, "distance", "matrix_binary_temp")
		if not os.path.exists(dir): os.mkdir(dir)

	def save(self):

		self.database_file = open(self.database_filename, "w")
		self.database_file.write("ID;Kmer_Spectrum_Filename\n")

		#for id, filename in self.items.items():
		for id in self.entries:
			filename = self.entries_infos[id]
			self.database_file.write(id + ";" + filename + "\n")

		self.database_file.close()

	def add_entry(self, id, relative_filename):
		self.entries.append(id)
		self.entries_infos[id] = relative_filename

		self.database_file = open(self.database_filename, "a")
		self.database_file.write(id + ";" + relative_filename + "\n")
		self.database_file.close()

	def contains_entry(self, id):
		return id in self.entries_infos

	#Called after a merge of k-mer spectrums, all merged dataset have to changed their path to the merged destination
	def change_entries(self, mergedDirs, mergeOutputRelativeDir):
		#new_id_filename = self.get_kmer_spectrum_dir_of_id(new_id, False)

		for id, filename in self.entries_infos.items():
			#filename_id = self.get_id_from_dir(filename)
			if filename in mergedDirs:
				self.entries_infos[id] = mergeOutputRelativeDir

	def get_default_kmer_spectrum_dir_of_id(self, id, abs_path):
		rel_path = os.path.join(self.kmer_spectrums_relative_dir, id)

		if abs_path:
			return os.path.join(self.dirname, rel_path)

		return rel_path

	def get_kmer_spectrum_dir_of_id(self, id, abs_path):
		rel_path = self.entries_infos[id]
		#rel_path = os.path.join(self.items[id], id)

		if abs_path:
			return os.path.join(self.dirname, rel_path)

		return rel_path

	def get_id_from_dir(self, dir):
		return os.path.basename(dir)


