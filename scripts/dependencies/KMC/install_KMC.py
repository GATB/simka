#Link to KMC binairies: http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download
#"""
#	Note:
#	For compilation on mac:
#	- kmc_api/kmer_defs.h : replace #include<ext/algorithm> by #include<algorithm>
#	- kmc_api/stdafx.h : replace #include<ext/algorithm> by #include<algorithm>
#	"""

#import os, shutil, sys
#os.chdir(os.path.split(os.path.realpath(__file__))[0])

import os
from sys import platform

os.chdir("scripts/simka2/bin")

if platform == "darwin":
	os.system("wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.mac.tar.gz")
	os.system("tar -zxf KMC3.mac.tar.gz")
	os.system("rm KMC3.mac.tar.gz")
else: #we expect linux
	os.system("wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz")
	os.system("tar -zxf KMC3.linux.tar.gz")
	os.system("rm KMC3.linux.tar.gz")
