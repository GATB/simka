
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include <string>
#include <iostream>
using namespace std;

int main (int argc, char* argv[])
{

	string command = "nohup ";

	for (int i = 1; i < argc; ++i) {
		//std::cout << argv[i] << std::endl;
		command += string(argv[i]) + " ";
	}

	//cout << command << endl;

	//cout << argc << " " << argv << endl;
	int ret=1;
	int nbTries = 0;
	while(ret != 0){
		ret = system(command.c_str());
		nanosleep((const struct timespec[]){{0, 10000000L}}, NULL);
		if(nbTries > 3) exit(1);
		nbTries += 1;
	}
}
