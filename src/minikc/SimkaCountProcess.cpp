
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
	int seg_fault_error = 11;
	int core_dump_error = 35584;
	int core_dump_error2 = 32256;
	int ret = core_dump_error;

	while(ret == core_dump_error || ret == core_dump_error2 || ret == seg_fault_error){
		ret = system(command.c_str());
		nanosleep((const struct timespec[]){{0, 10000000L}}, NULL);
	}
}
