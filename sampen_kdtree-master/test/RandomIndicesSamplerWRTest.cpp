#include <iostream>
#include <chrono>
#include <algorithm>

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <gsl/gsl_qrng.h>

#include "random_sampler.h"
#include "random_sampler.cpp"

//using namespace kdtree_mddc;
using std::cout;
using std::endl;
using std::string;
using std::vector;

int main(){
	vector<unsigned> random1 = GetSampleIndicesWR(SWR, 10000, 20, 5);
	cout << "size:" << random1.size() << endl;
	printf("random1结果如下：\n");
	for(int i = 0; i < random1.size(); i++)
	{
		cout << random1[i] << " ";
	}
	cout << endl;

}