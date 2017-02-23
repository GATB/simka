#Link to KMC binairies: http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download
"""
	Note:
	For compilation on mac:
	- kmc_api/kmer_defs.h : replace #include<ext/algorithm> by #include<algorithm>
	- kmc_api/stdafx.h : replace #include<ext/algorithm> by #include<algorithm>
	"""

import os, shutil, sys
os.chdir(os.path.split(os.path.realpath(__file__))[0])

#wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz
#tar -zxf KMC3.linux.tar.gz    -> dans le bin dir de simka!