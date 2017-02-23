

#Link to KMC binairies: http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download
"""
	Note:
		For compilation on mac:
		- kmc_api/kmer_defs.h : replace #include<ext/algorithm> by #include<algorithm>
		- kmc_api/stdafx.h : replace #include<ext/algorithm> by #include<algorithm>
"""

import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])

command = "wget https://github.com/refresh-bio/KMC/archive/v3.0.1.zip"